(* test_regression.sml -- OLS linear regression against hand-computed fits.

   Vector 1 (exact line y = 2x + 1): slope 2, intercept 1, r2 = 1.
   Vector 2: pts (1,1)(2,3)(3,2)(4,5)(5,4) with x-bar = y-bar = 3:
     Sxy = 8, Sxx = 10, Syy = 10  ->  slope 0.8, intercept 0.6, r2 = 0.64. *)

structure RegressionTests =
struct
  open Support
  structure S = Stats

  fun run () =
    let
      val () = Harness.section "linregress: exact line"
      val exact = [(0.0, 1.0), (1.0, 3.0), (2.0, 5.0), (3.0, 7.0), (4.0, 9.0)]
      val { slope, intercept, r2 } = S.linregress exact
      val () = checkApprox "exact slope = 2" (2.0, slope)
      val () = checkApprox "exact intercept = 1" (1.0, intercept)
      val () = checkApprox "exact r2 = 1" (1.0, r2)

      val () = Harness.section "linregress: known scatter"
      val pts = [(1.0, 1.0), (2.0, 3.0), (3.0, 2.0), (4.0, 5.0), (5.0, 4.0)]
      val r = S.linregress pts
      val () = checkApprox "scatter slope = 0.8" (0.8, #slope r)
      val () = checkApprox "scatter intercept = 0.6" (0.6, #intercept r)
      val () = checkApprox "scatter r2 = 0.64" (0.64, #r2 r)

      val () = Harness.section "linregress: negative slope"
      val dec = [(0.0, 10.0), (1.0, 8.0), (2.0, 6.0), (3.0, 4.0), (4.0, 2.0)]
      val rd = S.linregress dec
      val () = checkApprox "decreasing slope = -2" (~2.0, #slope rd)
      val () = checkApprox "decreasing intercept = 10" (10.0, #intercept rd)
      val () = checkApprox "decreasing r2 = 1" (1.0, #r2 rd)

      val () = Harness.section "linregress: edge cases"
      val () = Harness.checkRaises "single point raises"
                 (fn () => S.linregress [(1.0, 1.0)])
      val () = Harness.checkRaises "zero x-variance raises"
                 (fn () => S.linregress [(2.0, 1.0), (2.0, 3.0), (2.0, 5.0)])
    in
      ()
    end
end
