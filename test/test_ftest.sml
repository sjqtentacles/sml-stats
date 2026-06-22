(* test_ftest.sml -- F-distribution survival function and variance-ratio test.

   fSf(x,d1,d2) = P(F > x) via the regularized incomplete beta. Tabulated
   upper 5% critical values F_0.05(d1,d2) should map back to a survival of
   ~0.05 (|error| < 1e-3):
     F_0.05(5,5)   = 5.0503
     F_0.05(10,10) = 2.9782
     F_0.05(4,4)   = 6.3882
   fTest(a,b) statistic = variance a / variance b; dfn = |a|-1, dfd = |b|-1.
     a = [1,2,3,4,5]      -> sample variance 2.5
     b = [4,8,12,16,20]   -> sample variance 40.0
     statistic = 2.5/40 = 0.0625. *)

structure FTestTests =
struct
  open Support
  structure S = Stats

  val pTol = 1E~3

  fun run () =
    let
      val () = Harness.section "fSf survival vs F-tables (alpha = 0.05)"
      val () = checkApproxTol pTol "fSf(5.0503, 5, 5) ~ 0.05"
                 (0.05, S.fSf (5.0503, 5, 5))
      val () = checkApproxTol pTol "fSf(2.9782, 10, 10) ~ 0.05"
                 (0.05, S.fSf (2.9782, 10, 10))
      val () = checkApproxTol pTol "fSf(6.3882, 4, 4) ~ 0.05"
                 (0.05, S.fSf (6.3882, 4, 4))

      val () = Harness.section "fSf sanity"
      val () = Harness.check "fSf(1.0, 5, 5) in (0,1)"
                 (let val p = S.fSf (1.0, 5, 5) in p > 0.0 andalso p < 1.0 end)
      val () = Harness.check "fSf decreasing in x"
                 (S.fSf (1.0, 5, 7) > S.fSf (3.0, 5, 7))
      val () = checkApprox "fSf at 0 is 1" (1.0, S.fSf (0.0, 3, 3))

      val () = Harness.section "fTest variance ratio"
      val r0 = S.fTest ([1.0,2.0,3.0,4.0,5.0],[1.0,2.0,3.0,4.0,5.0])
      val () = checkApprox "statistic 1 for identical variances"
                 (1.0, #statistic r0)
      val () = Harness.checkInt "dfn = 4" (4, #dfn r0)
      val () = Harness.checkInt "dfd = 4" (4, #dfd r0)
      val r1 = S.fTest ([1.0,2.0,3.0,4.0,5.0],[4.0,8.0,12.0,16.0,20.0])
      val () = checkApprox "statistic = var a / var b = 2.5/40 = 0.0625"
                 (2.5 / 40.0, #statistic r1)
      val () = Harness.check "pValue in (0,1)"
                 (#pValue r1 > 0.0 andalso #pValue r1 < 1.0)
    in
      ()
    end
end
