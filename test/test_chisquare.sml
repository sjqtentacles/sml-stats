(* test_chisquare.sml -- chi-square goodness-of-fit test.

   statistic = Sum (o-e)^2/e; df = k-1; pValue = Q(df/2, statistic/2).

   Vectors:
     o = e = [10,10,10,10]      -> statistic 0, df 3, pValue 1.
     o = [16,18,16,14,12,12] (sum 88) vs expected 88/6 = 14.66667 each:
       Sum (o-e)^2 = 29.33333..., /14.66667 = 2.0 exactly; df = 5.
   Upper-tail vs tables (alpha = 0.05, |error| < 1e-3): the chi-square 5%
   critical values are 3.841 (df 1), 5.991 (df 2), 7.815 (df 3). Feeding a
   single-cell deviation of sqrt(target) (so statistic = target) over k cells
   (df = k-1) should give pValue ~ 0.05. *)

structure ChiSquareTests =
struct
  open Support
  structure S = Stats

  val pTol = 1E~3

  (* statistic = target with df = k-1, via one cell offset by sqrt target. *)
  fun chiP (target, k) =
    let
      val observed = (1.0 + Math.sqrt target) :: List.tabulate (k - 1, fn _ => 1.0)
      val expected = List.tabulate (k, fn _ => 1.0)
    in
      #pValue (S.chiSquareTest (observed, expected))
    end

  fun run () =
    let
      val () = Harness.section "chiSquareTest goodness-of-fit"
      val r0 = S.chiSquareTest ([10.0,10.0,10.0,10.0],[10.0,10.0,10.0,10.0])
      val () = checkApprox "statistic 0 when o = e" (0.0, #statistic r0)
      val () = Harness.checkInt "df = k-1 = 3" (3, #df r0)
      val () = checkApprox "pValue 1 when statistic 0" (1.0, #pValue r0)

      val e = 88.0 / 6.0
      val r1 = S.chiSquareTest ([16.0,18.0,16.0,14.0,12.0,12.0],
                                [e,e,e,e,e,e])
      val () = checkApproxTol 1E~6 "statistic = 2.0" (2.0, #statistic r1)
      val () = Harness.checkInt "df = 5" (5, #df r1)
      val () = Harness.check "pValue in (0,1)"
                 (#pValue r1 > 0.0 andalso #pValue r1 < 1.0)

      val () = Harness.section "chi-square upper tail vs tables (alpha = 0.05)"
      val () = checkApproxTol pTol "Q tail at 3.841, df 1 ~ 0.05"
                 (0.05, chiP (3.841, 2))
      val () = checkApproxTol pTol "Q tail at 5.991, df 2 ~ 0.05"
                 (0.05, chiP (5.991, 3))
      val () = checkApproxTol pTol "Q tail at 7.815, df 3 ~ 0.05"
                 (0.05, chiP (7.815, 4))
    in
      ()
    end
end
