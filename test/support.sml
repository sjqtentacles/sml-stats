(* support.sml -- shared helpers for the sml-stats tests.

   Every statistic is floating point, so comparisons go through an explicit
   epsilon (`approx`) rather than string or structural equality: `Real.toString`
   differs between MLton and Poly/ML, and the closed-form expectations only
   match up to rounding. A loose `eps` (1e-9) pins the algebraic identities;
   sampling-based checks use their own wider tolerances inline. *)

structure Support =
struct
  val eps = 1E~9

  fun approx (a, b) = Real.abs (a - b) <= eps

  (* approx with a caller-supplied tolerance, for Monte-Carlo style checks. *)
  fun approxTol tol (a, b) = Real.abs (a - b) <= tol

  fun checkApprox name (expected, actual) =
    Harness.check name (approx (expected, actual))

  fun checkApproxTol tol name (expected, actual) =
    Harness.check name (approxTol tol (expected, actual))

  (* Default generator state used throughout the suite. The default `Stats`
     structure is `StatsFn (SplitMix64)`, so its `rng` is `SplitMix64.state`. *)
  val seeded = SplitMix64.seed
end
