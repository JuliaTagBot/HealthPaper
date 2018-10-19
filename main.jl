function main(IG::InitialGuess,P::Params)

iter = 0
eflag = 0
max_iter = 1000
tol      = 1e-2
wgt      = 0.25

@everywhere K        =   IG.K
@everywhere Le        =   IG.L
@everywhere T       =  IG.T
@everywhere  τw       = IG.τw
while eflag==0

      @everywhere @unpack A,α,δ = P.F



      @everywhere w      = (1-α)*A*(K/Le)^α
      @everywhere r      = α*A*(K/Le)^(α-1)-δ

      @everywhere ss = P.G.ss_rate*w*Le

      println( "-------------------------")
      println( "-------------------------")
      @printf  " * w = %.3f\n" w
      @printf  " * r = %.3f\n" r
      @printf  " * τw = %.3f\n" τw
      println( "-------------------------")
      println( "-------------------------")

      #%%%%% Value function iteration
       vfi = VFI(IG,w, r, T, ss, τw,P)
      #%%%%% calculate invariant distribution


       mc= MC(vfi, w, r, ss, τw, T, P)

      # calculate output, given K = Khat
      Y = A*(mc.K^α)*(mc.Le^(1.0-α))

      err_K = abs(mc.K - K)
      err_Le = abs(mc.Le - Le)
      err_T = abs(mc.T - T)
      err_τw = abs(mc.τw - τw)

      iter = iter + 1

      println( "**************************")
      println( "iter = ", iter)
      println( "   K = ", K,    "   K' = ", mc.K)
      println( "   L = ", Le,    "   L' = ", mc.Le)
      println( "   T = ", T,    "   T' = ", mc.T)
      println( "tauw = ", τw, "   tauw' = ", mc.τw)
      println( " K/Y = ", mc.K/Y)
      println( " K/L = ", mc.K/mc.Le)
      println( " C/Y = ", mc.C/Y)
      println( " M/Y = ", mc.M/Y)
      println( "**************************")

      test1 = (err_K <= tol)
      test2 = (err_Le <= tol)
      test3 = (err_T <= tol)
      test4 = (err_τw <= tol)

      # for debugging purposes:
      #STOP

      if  test1 && test2 && test3 && test4
          eflag = 1
          println( "====================================")
          println( "CONVERGENCE")
          println( "====================================")
          #JLD.save("solution.jld","Khat",Khat, "Lhat",Lhat, "That",That, "tauwhat",tauwhat, "K_vec",K_vec,"L_vec",L_vec,"C_vec",C_vec,"F_vec",F_vec,"M_vec",M_vec,
          #"g_c",g_c,"v",v, "g_f",g_f, "g_n",g_n, "g_h",g_h, "g_an",g_an, "g_m",g_m)
      end

      if iter >= max_iter
          eflag = 2
          println( "====================================")
          println( "MAXIMUM NUMBER OF ITERATIONS REACHED")
          println( "====================================")
      end



      # update K,L,T,tauw
      K    = wgt*mc.K + (1.0-wgt)*K
      Le    = wgt*mc.Le + (1.0-wgt)*Le
      T    = wgt*mc.T + (1.0-wgt)*T
      τw = wgt*mc.τw + (1.0-wgt)*τw

end

      return mc, vfi

end

#end
#end
