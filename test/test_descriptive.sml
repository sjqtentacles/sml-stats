(* test_descriptive.sml -- closed-form descriptive-statistics vectors.

   The reference values are computed by hand for small fixed inputs:
     xs = [1,2,3,4,5]   mean 3, popVar 2, sampleVar 2.5
     ys = [1,2,3,4]     median 2.5 (even length) *)

structure DescriptiveTests =
struct
  open Support
  structure S = Stats

  val xs = [1.0, 2.0, 3.0, 4.0, 5.0]
  val ys = [1.0, 2.0, 3.0, 4.0]

  fun run () =
    let
      val () = Harness.section "sum / mean"
      val () = checkApprox "sum [1..5] = 15" (15.0, S.sum xs)
      val () = checkApprox "mean [1..5] = 3" (3.0, S.mean xs)
      val () = checkApprox "mean singleton" (7.0, S.mean [7.0])

      val () = Harness.section "variance / stddev"
      (* population variance of [1..5] = ((-2)^2+(-1)^2+0+1+4)/5 = 10/5 = 2 *)
      val () = checkApprox "variancePop [1..5] = 2" (2.0, S.variancePop xs)
      (* sample variance divides by n-1 = 4: 10/4 = 2.5 *)
      val () = checkApprox "variance (sample) [1..5] = 2.5" (2.5, S.variance xs)
      val () = checkApprox "stddevPop = sqrt 2"
                 (Math.sqrt 2.0, S.stddevPop xs)
      val () = checkApprox "stddev (sample) = sqrt 2.5"
                 (Math.sqrt 2.5, S.stddev xs)

      val () = Harness.section "median"
      val () = checkApprox "median odd [1..5] = 3" (3.0, S.median xs)
      val () = checkApprox "median even [1..4] = 2.5" (2.5, S.median ys)
      val () = checkApprox "median unsorted = 3"
                 (3.0, S.median [5.0, 1.0, 3.0, 2.0, 4.0])

      val () = Harness.section "quantile (type 7)"
      val () = checkApprox "q0   = min"  (1.0, S.quantile 0.0 xs)
      val () = checkApprox "q1   = max"  (5.0, S.quantile 1.0 xs)
      val () = checkApprox "q0.5 = median" (3.0, S.quantile 0.5 xs)
      (* h = (n-1)*q = 4*0.25 = 1.0 -> exactly xs[1] = 2.0 *)
      val () = checkApprox "q0.25 = 2.0" (2.0, S.quantile 0.25 xs)
      val () = checkApprox "q0.75 = 4.0" (4.0, S.quantile 0.75 xs)
      (* interpolated: [1,2,3,4], h = 3*0.5 = 1.5 -> between 2 and 3 = 2.5 *)
      val () = checkApprox "q0.5 even = 2.5" (2.5, S.quantile 0.5 ys)

      val () = Harness.section "minimum / maximum"
      val () = checkApprox "min [5,1,3,2,4] = 1"
                 (1.0, S.minimum [5.0, 1.0, 3.0, 2.0, 4.0])
      val () = checkApprox "max [5,1,3,2,4] = 5"
                 (5.0, S.maximum [5.0, 1.0, 3.0, 2.0, 4.0])

      val () = Harness.section "descriptive edge cases"
      val () = Harness.checkRaises "mean [] raises" (fn () => S.mean [])
      val () = Harness.checkRaises "median [] raises" (fn () => S.median [])
      val () = Harness.checkRaises "variance [] raises" (fn () => S.variance [])
    in
      ()
    end
end
