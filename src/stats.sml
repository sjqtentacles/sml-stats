(* stats.sml

   Implementation of `STATS` as a functor over a `sml-prng` RANDOM generator,
   plus a default `Stats` instantiated with SplitMix64.

   Everything is Basis-library `real`/`int` arithmetic threaded through pure
   helpers, so results are deterministic and identical under MLton and
   Poly/ML. The only non-elementary pieces are:
     - `erf` via the Abramowitz & Stegun 7.1.26 rational approximation, used by
       the normal cdf (|error| < 1.5e-7);
     - `gammaln` (Lanczos) + the regularized incomplete beta `betai`
       (Numerical Recipes continued fraction), used by the Student-t cdf and
       the t-test p-values. *)

functor StatsFn (R : RANDOM) :> STATS where type rng = R.state =
struct
  type rng = R.state

  exception Empty

  val pi = Math.pi

  (* ---- descriptive statistics ---- *)

  fun sum xs = List.foldl (op +) 0.0 xs

  fun count xs = real (List.length xs)

  fun mean xs =
    case xs of
      [] => raise Empty
    | _  => sum xs / count xs

  (* Sum of squared deviations from the mean. *)
  fun ss xs =
    let val m = mean xs
    in List.foldl (fn (x, acc) => acc + (x - m) * (x - m)) 0.0 xs end

  fun variancePop xs =
    case xs of [] => raise Empty | _ => ss xs / count xs

  fun variance xs =
    case xs of
      []  => raise Empty
    | [_] => raise Empty
    | _   => ss xs / (count xs - 1.0)

  fun stddevPop xs = Math.sqrt (variancePop xs)
  fun stddev xs = Math.sqrt (variance xs)

  (* Ascending merge sort on a real list (Basis-only, stable). *)
  fun sorted xs =
    let
      fun merge ([], ys) = ys
        | merge (xs, []) = xs
        | merge (x :: xs, y :: ys) =
            if x <= y then x :: merge (xs, y :: ys)
            else y :: merge (x :: xs, ys)
      fun split [] = ([], [])
        | split [a] = ([a], [])
        | split (a :: b :: rest) =
            let val (l, r) = split rest in (a :: l, b :: r) end
      fun msort [] = []
        | msort [a] = [a]
        | msort ys =
            let val (l, r) = split ys
            in merge (msort l, msort r) end
    in
      msort xs
    end

  fun minimum xs =
    case xs of [] => raise Empty | y :: ys => List.foldl Real.min y ys
  fun maximum xs =
    case xs of [] => raise Empty | y :: ys => List.foldl Real.max y ys

  (* Linear-interpolation quantile (R/NumPy "type 7"). *)
  fun quantile q xs =
    case xs of
      [] => raise Empty
    | _  =>
        let
          val arr = Array.fromList (sorted xs)
          val n = Array.length arr
          val qq = if q < 0.0 then 0.0 else if q > 1.0 then 1.0 else q
          val h = real (n - 1) * qq
          val lo = Real.floor h
          val frac = h - real lo
          val xl = Array.sub (arr, lo)
        in
          if lo + 1 < n
          then xl + frac * (Array.sub (arr, lo + 1) - xl)
          else xl
        end

  fun median xs = quantile 0.5 xs

  (* ---- special functions ---- *)

  (* erf via Abramowitz & Stegun 7.1.26 (|error| < 1.5e-7). *)
  fun erf x =
    let
      val sign = if x < 0.0 then ~1.0 else 1.0
      val ax = Real.abs x
      val t = 1.0 / (1.0 + 0.3275911 * ax)
      val poly =
        ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t
          - 0.284496736) * t + 0.254829592) * t
      val y = 1.0 - poly * Math.exp (~(ax * ax))
    in
      sign * y
    end

  (* Lanczos log-gamma, g = 7, n = 9 coefficients. *)
  local
    val g = 7.0
    val c =
      [ 0.99999999999980993,
        676.5203681218851,
        ~1259.1392167224028,
        771.32342877765313,
        ~176.61502916214059,
        12.507343278686905,
        ~0.13857109526572012,
        9.9843695780195716E~6,
        1.5056327351493116E~7 ]
  in
    fun gammaln z0 =
      let
        val z = z0 - 1.0
        val (a, _) =
          List.foldl
            (fn (ci, (acc, i)) =>
               if i = 0 then (ci, 1)
               else (acc + ci / (z + real i), i + 1))
            (0.0, 0) c
        val tt = z + g + 0.5
      in
        0.5 * Math.ln (2.0 * pi) + (z + 0.5) * Math.ln tt - tt + Math.ln a
      end
  end

  (* Continued fraction for the incomplete beta (Numerical Recipes betacf). *)
  fun betacf (a, b, x) =
    let
      val maxIt = 200
      val epsCf = 3.0E~12
      val fpmin = 1.0E~300
      val qab = a + b
      val qap = a + 1.0
      val qam = a - 1.0
      fun guard d = if Real.abs d < fpmin then fpmin else d
      fun loop (m, c, d, h) =
        if m > maxIt then h
        else
          let
            val rm = real m
            (* even step *)
            val aa1 = rm * (b - rm) * x / ((qam + 2.0 * rm) * (a + 2.0 * rm))
            val d1 = guard (1.0 + aa1 * d)
            val c1 = guard (1.0 + aa1 / c)
            val d1' = 1.0 / d1
            val h1 = h * d1' * c1
            (* odd step *)
            val aa2 =
              ~((a + rm) * (qab + rm) * x
                / ((a + 2.0 * rm) * (qap + 2.0 * rm)))
            val d2 = guard (1.0 + aa2 * d1')
            val c2 = guard (1.0 + aa2 / c1)
            val d2' = 1.0 / d2
            val del = d2' * c2
            val h2 = h1 * del
          in
            if Real.abs (del - 1.0) < epsCf then h2
            else loop (m + 1, c2, d2', h2)
          end
      val d0 = guard (1.0 - qab * x / qap)
      val d0' = 1.0 / d0
    in
      loop (1, 1.0, d0', d0')
    end

  (* Regularized incomplete beta I_x(a,b). *)
  fun betai (a, b, x) =
    if x <= 0.0 then 0.0
    else if x >= 1.0 then 1.0
    else
      let
        val bt =
          Math.exp
            (gammaln (a + b) - gammaln a - gammaln b
             + a * Math.ln x + b * Math.ln (1.0 - x))
      in
        if x < (a + 1.0) / (a + b + 2.0)
        then bt * betacf (a, b, x) / a
        else 1.0 - bt * betacf (b, a, 1.0 - x) / b
      end

  (* ---- distributions ---- *)

  structure Normal =
  struct
    type param = { mu : real, sigma : real }

    fun pdf { mu, sigma } x =
      let val z = (x - mu) / sigma
      in Math.exp (~0.5 * z * z) / (sigma * Math.sqrt (2.0 * pi)) end

    fun cdf { mu, sigma } x =
      0.5 * (1.0 + erf ((x - mu) / (sigma * Math.sqrt 2.0)))

    (* Box-Muller: two uniforms -> one standard normal. *)
    fun sample { mu, sigma } s =
      let
        val (u1, s1) = R.real01 s
        val (u2, s2) = R.real01 s1
        (* guard ln 0 at the open end of [0,1) *)
        val u1' = if u1 <= 0.0 then 1.0E~300 else u1
        val z = Math.sqrt (~2.0 * Math.ln u1') * Math.cos (2.0 * pi * u2)
      in
        (mu + sigma * z, s2)
      end
  end

  structure Binomial =
  struct
    type param = { n : int, p : real }

    (* C(n,k) computed multiplicatively in reals (n is small in practice). *)
    fun choose (n, k) =
      if k < 0 orelse k > n then 0.0
      else
        let
          val k = Int.min (k, n - k)
          fun loop (i, acc) =
            if i > k then acc
            else loop (i + 1, acc * real (n - k + i) / real i)
        in
          loop (1, 1.0)
        end

    fun pdf { n, p } k =
      if k < 0 orelse k > n then 0.0
      else choose (n, k) * Math.pow (p, real k)
           * Math.pow (1.0 - p, real (n - k))

    fun cdf (par as { n, ... }) k =
      if k < 0 then 0.0
      else
        let
          val kk = Int.min (k, n)
          fun loop (i, acc) =
            if i > kk then acc else loop (i + 1, acc + pdf par i)
        in
          loop (0, 0.0)
        end

    fun sample { n, p } s =
      let
        fun loop (0, acc, s) = (acc, s)
          | loop (i, acc, s) =
              let val (u, s') = R.real01 s
              in loop (i - 1, if u < p then acc + 1 else acc, s') end
      in
        loop (n, 0, s)
      end
  end

  structure Poisson =
  struct
    type param = { lambda : real }

    fun pdf { lambda } k =
      if k < 0 then 0.0
      else Math.exp (~lambda + real k * Math.ln lambda - gammaln (real k + 1.0))

    fun cdf (par as { lambda }) k =
      if k < 0 then 0.0
      else
        let
          fun loop (i, acc) =
            if i > k then acc else loop (i + 1, acc + pdf par i)
        in
          loop (0, 0.0)
        end

    (* Knuth's multiplicative algorithm. *)
    fun sample { lambda } s =
      let
        val lcap = Math.exp (~lambda)
        fun loop (k, p, s) =
          let
            val (u, s') = R.real01 s
            val p' = p * u
          in
            if p' <= lcap then (k, s') else loop (k + 1, p', s')
          end
      in
        loop (0, 1.0, s)
      end
  end

  (* ---- ordinary least-squares regression ---- *)

  fun linregress pts =
    let
      val n = real (List.length pts)
      val () = if n < 2.0 then raise Empty else ()
      val (sx, sy) =
        List.foldl (fn ((x, y), (ax, ay)) => (ax + x, ay + y)) (0.0, 0.0) pts
      val mx = sx / n
      val my = sy / n
      val (sxx, syy, sxy) =
        List.foldl
          (fn ((x, y), (axx, ayy, axy)) =>
             let val dx = x - mx val dy = y - my
             in (axx + dx * dx, ayy + dy * dy, axy + dx * dy) end)
          (0.0, 0.0, 0.0) pts
      val () = if sxx <= 0.0 then raise Empty else ()
      val slope = sxy / sxx
      val intercept = my - slope * mx
      val r2 = if syy <= 0.0 then 1.0 else (sxy * sxy) / (sxx * syy)
    in
      { slope = slope, intercept = intercept, r2 = r2 }
    end

  (* ---- Student's t ---- *)

  (* CDF via the incomplete beta: with x = df/(df+t^2),
       ib = I_x(df/2, 1/2) = P(|T| >= |t|)  (the two-tailed tail mass).
     Then cdf(t) = 1 - ib/2 for t > 0, and ib/2 for t <= 0. *)
  fun studentTCdf { t, df } =
    let
      val x = df / (df + t * t)
      val ib = betai (df / 2.0, 0.5, x)
    in
      if t > 0.0 then 1.0 - 0.5 * ib else 0.5 * ib
    end

  (* Two-tailed p-value for a t statistic: I_x(df/2, 1/2), x = df/(df+t^2). *)
  fun twoTailedP (t, df) = betai (df / 2.0, 0.5, df / (df + t * t))

  fun tTestOne { data, mu0 } =
    let
      val n = real (List.length data)
      val () = if n < 2.0 then raise Empty else ()
      val m = mean data
      val se = stddev data / Math.sqrt n
      val t = (m - mu0) / se
      val df = n - 1.0
    in
      { t = t, df = df, pValue = twoTailedP (t, df) }
    end

  fun tTestTwo { a, b } =
    let
      val n1 = real (List.length a)
      val n2 = real (List.length b)
      val () = if n1 < 2.0 orelse n2 < 2.0 then raise Empty else ()
      val m1 = mean a
      val m2 = mean b
      val v1 = variance a
      val v2 = variance b
      val df = n1 + n2 - 2.0
      val sp2 = ((n1 - 1.0) * v1 + (n2 - 1.0) * v2) / df
      val se = Math.sqrt (sp2 * (1.0 / n1 + 1.0 / n2))
      val t = (m1 - m2) / se
    in
      { t = t, df = df, pValue = twoTailedP (t, df) }
    end
end

(* Default instantiation: deterministic SplitMix64 sampling. *)
structure Stats = StatsFn (SplitMix64)
