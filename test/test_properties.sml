(* test_properties.sml -- property-based tests (sml-check) for sml-stats.

   These properties check algebraic invariants that must hold for ANY input,
   as opposed to the closed-form checks in the other suites which pin
   specific hand-computed values:

     - the mean of a constant list is that constant;
     - the (sample) variance and stddev of a constant list are 0;
     - the mean is invariant under reordering the input (it is sum-based,
       so summation order shouldn't matter up to floating-point rounding);
     - the median always lies between the minimum and the maximum;
     - a non-constant sample is perfectly correlated with itself
       (correlation(x,x) = 1).

   All comparisons go through an explicit epsilon per Support's convention;
   `Real.toString` is never used for display, only the fixed-decimal
   `showReal` below (its output only actually reaches stdout on a shrunk
   counterexample). *)

structure PropertyTests =
struct
  open Support
  structure S = Stats

  (* Fixed-decimal real formatting for shrunk-counterexample display; never
     Real.toString (its rendering differs between MLton and Poly/ML).
     Negative zero is normalized to positive. *)
  fun showReal x =
    let val x = if Real.== (x, 0.0) then 0.0 else x
    in Real.fmt (StringCvt.FIX (SOME 6)) x end

  fun showRealList xs = "[" ^ String.concatWith "," (List.map showReal xs) ^ "]"

  fun run () =
    let
      val () = Harness.section "properties: descriptive statistics (sml-check)"

      val genVal = Check.realRange (~1000.0, 1000.0)
      val genLen = Check.choose (1, 30)
      val genLen2 = Check.choose (2, 30)

      (* ---- mean of a constant list equals that constant ---- *)
      val genConst1 =
        Check.bind genVal (fn c =>
          Check.map (fn n => (c, n)) genLen)
      fun showConst (c, n) = showReal c ^ " x" ^ Int.toString n

      val () =
        Harness.check "prop: mean of a constant list equals that constant"
          (case Check.quickCheck
                  (Check.forAll genConst1 showConst
                     (fn (c, n) =>
                        let val xs = List.tabulate (n, fn _ => c)
                        in approx (c, S.mean xs) end)) of
               Check.Passed _ => true
             | Check.Failed _ => false)

      (* ---- variance/stddev of a constant list are 0 (n >= 2, so the
         sample variance's n-1 divisor is well-defined) ---- *)
      val genConst2 =
        Check.bind genVal (fn c =>
          Check.map (fn n => (c, n)) genLen2)

      val () =
        Harness.check "prop: variance and stddev of a constant list are 0"
          (case Check.quickCheck
                  (Check.forAll genConst2 showConst
                     (fn (c, n) =>
                        let val xs = List.tabulate (n, fn _ => c)
                        in approx (0.0, S.variance xs)
                           andalso approx (0.0, S.stddev xs)
                           andalso approx (0.0, S.variancePop xs)
                           andalso approx (0.0, S.stddevPop xs)
                        end)) of
               Check.Passed _ => true
             | Check.Failed _ => false)

      (* ---- mean is invariant under reversal (sum-based, so order
         shouldn't matter beyond floating-point rounding) ---- *)
      val genList = Check.nonEmptyListOf genVal

      val () =
        Harness.check "prop: mean is invariant under list reversal"
          (case Check.quickCheck
                  (Check.forAll genList showRealList
                     (fn xs => approx (S.mean xs, S.mean (List.rev xs)))) of
               Check.Passed _ => true
             | Check.Failed _ => false)

      (* ---- median always lies between the minimum and the maximum ---- *)
      val () =
        Harness.check "prop: median lies between minimum and maximum"
          (case Check.quickCheck
                  (Check.forAll genList showRealList
                     (fn xs =>
                        let
                          val m  = S.median xs
                          val lo = S.minimum xs
                          val hi = S.maximum xs
                        in
                          m >= lo - eps andalso m <= hi + eps
                        end)) of
               Check.Passed _ => true
             | Check.Failed _ => false)

      (* ---- a non-constant sample is perfectly correlated with itself.
         Continuous reals are non-constant with probability 1, so the
         filter below essentially never resamples. ---- *)
      val genNonConstList =
        Check.filter
          (fn xs => not (List.all (fn y => Real.abs (y - List.hd xs) < 1E~9) xs))
          (Check.bind genLen2 (fn n => Check.listOfLen n genVal))

      val () =
        Harness.check "prop: correlation(x,x) = 1 for a non-constant sample"
          (case Check.quickCheck
                  (Check.forAll genNonConstList showRealList
                     (fn xs => approx (1.0, S.correlation (xs, xs)))) of
               Check.Passed _ => true
             | Check.Failed _ => false)
    in
      ()
    end
end
