(* test_ttest.sml -- Student-t CDF and t-tests against published values.

   studentTCdf vectors:
     cdf(t=0,  df=5)  = 0.5            (symmetry)
     cdf(t=1,  df=1)  = 0.75           (Cauchy: 0.5 + atan(1)/pi)
     cdf(t=-1, df=1)  = 0.25
   t-table vectors:
     df=10, |t|=2.228 -> two-tailed p ~ 0.05
   Two-sample pooled (clean by construction):
     a=[1..5], b=[6..10]: m1=3, m2=8, v1=v2=2.5, sp2=2.5,
     se = sqrt(2.5*(2/5)) = 1, t = (3-8)/1 = -5, df = 8. *)

structure TTestTests =
struct
  open Support
  structure S = Stats

  val tTol = 1E~5

  fun run () =
    let
      val () = Harness.section "studentTCdf closed forms"
      val () = checkApproxTol tTol "cdf(0, df=5) = 0.5"
                 (0.5, S.studentTCdf { t = 0.0, df = 5.0 })
      val () = checkApproxTol tTol "cdf(1, df=1) = 0.75 (Cauchy)"
                 (0.75, S.studentTCdf { t = 1.0, df = 1.0 })
      val () = checkApproxTol tTol "cdf(-1, df=1) = 0.25"
                 (0.25, S.studentTCdf { t = ~1.0, df = 1.0 })
      (* Cauchy closed form at several points: F(t) = 0.5 + atan(t)/pi *)
      val () = checkApproxTol tTol "cdf(2, df=1) = 0.5+atan2/pi"
                 (0.5 + Math.atan 2.0 / Math.pi,
                  S.studentTCdf { t = 2.0, df = 1.0 })
      (* monotone, in (0,1) *)
      val () = Harness.check "cdf increasing in t"
                 (S.studentTCdf { t = 1.0, df = 8.0 }
                  < S.studentTCdf { t = 2.0, df = 8.0 })

      val () = Harness.section "one-sample t-test"
      (* data = [1..5], mu0 = 3 = mean -> t = 0, p = 1 *)
      val r0 = S.tTestOne { data = [1.0,2.0,3.0,4.0,5.0], mu0 = 3.0 }
      val () = checkApproxTol tTol "t = 0 when mu0 = mean" (0.0, #t r0)
      val () = checkApproxTol tTol "df = n-1 = 4" (4.0, #df r0)
      val () = checkApproxTol tTol "p = 1 when t = 0" (1.0, #pValue r0)
      (* data = [1..5], mu0 = 0: mean 3, s = sqrt 2.5, se = s/sqrt 5,
         t = 3 / (sqrt2.5/sqrt5) = 3 / sqrt(0.5) = 4.2426407 *)
      val r1 = S.tTestOne { data = [1.0,2.0,3.0,4.0,5.0], mu0 = 0.0 }
      val () = checkApproxTol tTol "one-sample t = 4.2426407"
                 (3.0 / Math.sqrt 0.5, #t r1)
      val () = Harness.check "one-sample p < 0.05" (#pValue r1 < 0.05)

      val () = Harness.section "t-table p-value (df=10, t=2.228)"
      (* critical t for two-tailed alpha=0.05, df=10 is 2.228; p ~ 0.05 *)
      val () = checkApproxTol 1E~3 "p(|t|=2.228, df=10) ~ 0.05"
                 (0.05, S.studentTCdf { t = ~2.228, df = 10.0 } * 2.0)

      val () = Harness.section "two-sample pooled t-test"
      val r2 = S.tTestTwo { a = [1.0,2.0,3.0,4.0,5.0],
                            b = [6.0,7.0,8.0,9.0,10.0] }
      val () = checkApproxTol tTol "two-sample t = -5" (~5.0, #t r2)
      val () = checkApproxTol tTol "two-sample df = 8" (8.0, #df r2)
      (* t=-5, df=8 -> two-tailed p ~ 0.00105 *)
      val () = Harness.check "two-sample p very small" (#pValue r2 < 0.002)
      val () = checkApproxTol 1E~4 "two-sample p ~ 0.00105"
                 (0.00105, #pValue r2)
      (* identical samples -> t = 0, p = 1 *)
      val r3 = S.tTestTwo { a = [1.0,2.0,3.0,4.0,5.0],
                            b = [1.0,2.0,3.0,4.0,5.0] }
      val () = checkApproxTol tTol "identical samples t = 0" (0.0, #t r3)
      val () = checkApproxTol tTol "identical samples p = 1" (1.0, #pValue r3)
    in
      ()
    end
end
