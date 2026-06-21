# sml-stats

Statistics in pure Standard ML — descriptive moments, three probability
distributions (Normal, Binomial, Poisson) with `pdf`/`cdf`/`sample`, ordinary
least-squares linear regression, and Student's t-test — built on
[`sml-prng`](https://github.com/sjqtentacles/sml-prng) for deterministic,
seedable sampling. No FFI, no external dependencies, and **deterministic**,
byte-identically under both [MLton](http://mlton.org/) and
[Poly/ML](https://www.polyml.org/).

## Status

- 78 assertions, green on MLton and Poly/ML.
- Basis-library + vendored `sml-prng` only; deterministic across compilers.
- Vendors `sml-prng` (Layout B), so the repo builds standalone.

## Install

With [`smlpkg`](https://github.com/diku-dk/smlpkg):

```
smlpkg add github.com/sjqtentacles/sml-stats
smlpkg sync
```

Include the MLB from your own (it pulls in the vendored `sml-prng`):

```
local
  $(SML_LIB)/basis/basis.mlb
  lib/github.com/sjqtentacles/sml-stats/... (via smlpkg)
in
  ...
end
```

This brings `structure Stats` (and the vendored generators) into scope.

## Quick start

```sml
(* descriptive statistics *)
val xs = [1.0, 2.0, 3.0, 4.0, 5.0]
val m  = Stats.mean xs              (* 3.0   *)
val sd = Stats.stddev xs            (* sqrt 2.5, sample stddev *)
val md = Stats.median xs            (* 3.0   *)
val q1 = Stats.quantile 0.25 xs     (* 2.0   *)

(* distributions: pdf / cdf *)
val p  = Stats.Normal.cdf { mu = 0.0, sigma = 1.0 } 1.96   (* ~0.975 *)
val pk = Stats.Binomial.pdf { n = 10, p = 0.5 } 5          (* 0.24609375 *)

(* seedable, pure sampling (thread the generator state) *)
val (z, s1) = Stats.Normal.sample { mu = 0.0, sigma = 1.0 } (SplitMix64.seed 0w42)

(* ordinary least-squares regression *)
val { slope, intercept, r2 } =
  Stats.linregress [(0.0,1.0),(1.0,3.0),(2.0,5.0),(3.0,7.0)]  (* 2.0, 1.0, 1.0 *)

(* Student's t-test *)
val { t, df, pValue } =
  Stats.tTestTwo { a = [1.0,2.0,3.0,4.0,5.0], b = [6.0,7.0,8.0,9.0,10.0] }
(* t = -5.0, df = 8.0 *)
```

## API (`signature STATS`)

```sml
(* descriptive (input lists must be non-empty; raise Empty otherwise) *)
val mean        : real list -> real
val variance    : real list -> real    (* sample,     /(n-1) *)
val variancePop : real list -> real    (* population, /n     *)
val stddev      : real list -> real
val stddevPop   : real list -> real
val median      : real list -> real
val quantile    : real -> real list -> real   (* type-7 linear interp *)
val minimum     : real list -> real
val maximum     : real list -> real

(* distributions: each has pdf, cdf, and a pure `sample` (rng -> value * rng) *)
structure Normal   : sig type param = { mu : real, sigma : real } ... end
structure Binomial : sig type param = { n : int, p : real }       ... end
structure Poisson  : sig type param = { lambda : real }           ... end

(* regression and inference *)
val linregress  : (real * real) list
                  -> { slope : real, intercept : real, r2 : real }
val studentTCdf : { t : real, df : real } -> real
val tTestOne    : { data : real list, mu0 : real }
                  -> { t : real, df : real, pValue : real }
val tTestTwo    : { a : real list, b : real list }
                  -> { t : real, df : real, pValue : real }
```

The library is a functor `StatsFn (R : RANDOM) :> STATS where type rng = R.state`
over a [`sml-prng`](https://github.com/sjqtentacles/sml-prng) generator. The
default `structure Stats = StatsFn (SplitMix64)`, so its `rng` is
`SplitMix64.state`; instantiate `StatsFn` yourself to sample from
`Xoshiro256ss` or `Pcg32` instead.

### Conventions

- `variance`/`stddev` are the **sample** (Bessel-corrected, `/(n-1)`)
  statistics; `variancePop`/`stddevPop` divide by `n`.
- `quantile q xs` uses the R/NumPy **type-7** linear-interpolation rule
  (`h = (n-1)·q`), so `quantile 0.5 = median`.
- `Normal.cdf` uses the Abramowitz & Stegun 7.1.26 `erf` approximation
  (|error| < 1.5e-7).
- `studentTCdf` and the t-test p-values use the regularized incomplete beta
  function (`gammaln` + a Numerical-Recipes continued fraction), matching
  published t-tables; the reported p-values are **two-tailed**.
- Sampling is **pure and seedable**: `sample` takes a generator state and
  returns the draw plus the successor state. The same seed yields the same
  draws on every run, machine, and compiler.

## Build & test

```
make test        # MLton
make test-poly   # Poly/ML
make all-tests   # both
make example     # build + run examples/demo.sml
make clean
```

Both compilers run the same strict-TDD suite, seeded with closed-form vectors:
hand-computed moments and quantiles, distribution `pdf`/`cdf` values
(`N(0,1)` cdf at 0/1/1.96, `Binom(10,½)` at `k = 5`, `Poisson(3)`), exact and
noisy regression fits, and t-statistics against known values (a two-sample
`t = -5, df = 8` by construction; the `df = 10, |t| = 2.228 → p ≈ 0.05`
t-table entry; the Cauchy `df = 1` closed form). Samplers are checked for
reproducibility, support bounds, and the law of large numbers.

## Example

`make example` draws 1000 `Normal(50, 10)` samples from a fixed seed and prints
a reproducible summary, an ASCII histogram, a regression fit, and two t-tests
(output is byte-identical under MLton and Poly/ML):

```
=== sml-stats demo ============================================

Drew 1000 samples from Normal(mu=50, sigma=10)
  with SplitMix64 seed 0x20260621 (deterministic).

Descriptive statistics
  mean      = 49.7818
  stddev    = 9.9300
  variance  = 98.6046
  min       = 18.0686
  q25       = 43.1503
  median    = 50.0329
  q75       = 56.6631
  max       = 80.4028

Histogram (Normal samples)
     18.1 | # 4
     20.7 | ## 5
     23.3 | #### 9
     25.9 | ## 5
     28.5 | ####### 15
     31.1 | ####### 15
     33.7 | ########### 23
     36.2 | ######################## 50
     38.8 | ############################### 63
     41.4 | ############################################ 89
     44.0 | ################################################ 98
     46.6 | ############################################## 94
     49.2 | ################################################# 99
     51.8 | ################################################## 101
     54.4 | ############################################## 93
     57.0 | ####################################### 80
     59.6 | ######################## 50
     62.2 | ######################### 51
     64.8 | ########## 21
     67.4 | ######### 20
     70.0 | ### 7
     72.6 | # 3
     75.2 | # 4
     77.8 |  1

Linear regression  (true model: y = 3x + 2 + N(0,2) noise)
  slope     = 3.0243   (true 3.0)
  intercept = 1.2239   (true 2.0)
  r^2       = 0.9899

One-sample t-test  (group A vs mu0 = 5.0)
  t = 0.8799   df = 6.0   p = 0.4128
Two-sample t-test  (group A vs group B, pooled)
  t = ~6.1619   df = 12.0   p = 0.000049

===============================================================
```

The PNG histogram/regression chart via `sml-plot` is deferred until that
library lands; the ASCII histogram above is the committed text asset.

### Poly/ML note

CI builds Poly/ML 5.9.1 from source rather than using the Ubuntu package
(Poly/ML 5.7.1), whose X86 code generator crashes (`asGenReg raised while
compiling`) on heavy real-arithmetic code. See `.github/workflows/ci.yml`.

## License

MIT — see [LICENSE](LICENSE).
