(* stats.sig

   A small, pure statistics toolkit: descriptive moments, three probability
   distributions (Normal, Binomial, Poisson) with pdf/cdf/sample, ordinary
   least-squares linear regression, and Student's t-test.

   Everything is built on Basis-library `real`/`int` and the vendored
   `sml-prng` generators, so it is deterministic and behaves identically under
   MLton and Poly/ML. Sampling is pure: every `sample` takes a generator state
   and returns the draw together with the successor state, so callers thread
   the state exactly as they would with `sml-prng`.

   The abstract type `rng` is the random-generator state. The library is a
   functor over a `RANDOM` generator (see `StatsFn`); the default `Stats`
   structure is instantiated with `SplitMix64`, so `rng = SplitMix64.state`.

   Numerical conventions:
     - `variance`/`stddev` are the *sample* statistics (Bessel-corrected,
       divide by n-1); `variancePop`/`stddevPop` are the population forms
       (divide by n).
     - `quantile q xs` uses the linear-interpolation method (R/NumPy "type 7":
       h = (n-1)*q), so `quantile 0.5` agrees with `median`.
     - `Normal.cdf` uses the Abramowitz & Stegun 7.1.26 `erf` approximation
       (|error| < 1.5e-7).
     - `studentTCdf`/the t-test p-values use the regularized incomplete beta
       function, so they match published t-tables. *)

signature STATS =
sig
  (* Generator state of the underlying `sml-prng` instance. *)
  type rng

  (* Raised by the descriptive statistics on an empty input. *)
  exception Empty

  (* ---- descriptive statistics (input lists must be non-empty) ---- *)
  val sum         : real list -> real
  val mean        : real list -> real
  val variance    : real list -> real   (* sample,     /(n-1) *)
  val variancePop : real list -> real   (* population, /n     *)
  val stddev      : real list -> real   (* sqrt variance      *)
  val stddevPop   : real list -> real   (* sqrt variancePop   *)
  val median      : real list -> real
  (* `quantile q xs`, q in [0,1], linear interpolation (type 7). *)
  val quantile    : real -> real list -> real
  val minimum     : real list -> real
  val maximum     : real list -> real

  (* ---- distributions ---- *)

  structure Normal :
  sig
    type param = { mu : real, sigma : real }     (* sigma > 0 *)
    val pdf    : param -> real -> real
    val cdf    : param -> real -> real
    (* Box-Muller draw. *)
    val sample : param -> rng -> real * rng
  end

  structure Binomial :
  sig
    type param = { n : int, p : real }           (* n >= 0, 0 <= p <= 1 *)
    val pdf    : param -> int -> real             (* P(X = k)            *)
    val cdf    : param -> int -> real             (* P(X <= k)           *)
    (* Sum of n Bernoulli(p) trials. *)
    val sample : param -> rng -> int * rng
  end

  structure Poisson :
  sig
    type param = { lambda : real }               (* lambda > 0 *)
    val pdf    : param -> int -> real             (* P(X = k)   *)
    val cdf    : param -> int -> real             (* P(X <= k)  *)
    (* Knuth's multiplicative algorithm. *)
    val sample : param -> rng -> int * rng
  end

  (* ---- ordinary least-squares regression ---- *)
  (* `linregress pts` fits y = slope*x + intercept; `r2` is the coefficient of
     determination. Requires at least two points with non-zero x-variance. *)
  val linregress :
    (real * real) list -> { slope : real, intercept : real, r2 : real }

  (* ---- correlation ---- *)

  (* Pearson product-moment correlation coefficient r of two equal-length,
     non-constant samples (xs, ys). Raises `Empty` if the lists differ in
     length, have fewer than two points, or either has zero variance.
     `correlation` is an alias for `pearson`. *)
  val pearson     : real list * real list -> real
  val correlation : real list * real list -> real

  (* Spearman rank correlation coefficient: Pearson's r computed on the
     fractional ranks of each sample, where tied values share their average
     rank. Same domain restrictions as `pearson`. *)
  val spearman    : real list * real list -> real

  (* ---- chi-square goodness-of-fit ---- *)

  (* `chiSquareTest (observed, expected)` for two equal-length lists of
     category counts (expected entries must be > 0). Returns
       statistic = Sum_i (o_i - e_i)^2 / e_i,
       df        = length - 1,
       pValue    = Q(df/2, statistic/2)   (regularized upper incomplete gamma).
     Raises `Empty` on length mismatch or fewer than two categories. *)
  val chiSquareTest :
    real list * real list -> { statistic : real, df : int, pValue : real }

  (* ---- F distribution ---- *)

  (* Survival function of the F distribution: `fSf (x, dfn, dfd)` = P(F > x)
     for an F variate with `dfn` numerator and `dfd` denominator degrees of
     freedom, via the regularized incomplete beta. `x` >= 0, dfn,dfd >= 1. *)
  val fSf : real * int * int -> real

  (* `fTest (a, b)` compares the two samples' variances with an F-test:
       statistic = variance a / variance b,
       dfn       = length a - 1,  dfd = length b - 1,
       pValue    = fSf (statistic, dfn, dfd)   (upper-tail).
     Both samples need at least two points with non-zero variance in `b`. *)
  val fTest :
    real list * real list ->
    { statistic : real, dfn : int, dfd : int, pValue : real }

  (* ---- Student's t ---- *)

  (* CDF of the Student-t distribution with `df` degrees of freedom. *)
  val studentTCdf : { t : real, df : real } -> real

  (* One-sample t-test of `data` against the hypothesized mean `mu0`. *)
  val tTestOne :
    { data : real list, mu0 : real } ->
    { t : real, df : real, pValue : real }

  (* Two-sample, equal-variance (pooled) Student's t-test. *)
  val tTestTwo :
    { a : real list, b : real list } ->
    { t : real, df : real, pValue : real }
end
