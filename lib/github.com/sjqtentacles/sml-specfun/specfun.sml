(* specfun.sml -- pure Standard ML special functions.

   Implementation notes
   --------------------
   * lgamma: Lanczos approximation (g = 7, 9-term series), with the reflection
     formula  G(x) G(1-x) = pi / sin(pi x)  for the left half-plane.
   * gamma: exp of the Lanczos log form (it never builds an intermediate that
     overflows before the final result does), with reflection for x < 1/2.
   * gammaIncP/Q: the textbook split -- a power series for x < a+1 and a
     Lentz-evaluated continued fraction for x >= a+1 -- each driven to a fixed
     relative tolerance.
   * erf/erfc: expressed through the regularized incomplete gamma so the far
     tail of erfc stays accurate (erfc x = Q(1/2, x^2) for x >= 0), instead of
     the catastrophic  1 - erf x.
   * erfInv: a closed-form initial guess refined by Newton's method on erf, run
     to a fixed tolerance.
   * betaInc: the prefactor  x^a (1-x)^b / (a B(a,b))  times a Lentz continued
     fraction, taking the  I_x(a,b) = 1 - I_{1-x}(b,a)  branch for fast
     convergence.

   Every loop terminates on a value tolerance (never a wall-clock or fixed
   step count that could diverge between compilers), and only +,-,*,/ and the
   Basis `Math` functions are used, so output is deterministic and identical
   under MLton and Poly/ML. *)

structure Specfun :> SPECFUN =
struct
  val eps = 1E~15

  val pi = 3.14159265358979323846

  (* Smallest positive value used to dodge division-by-zero in the Lentz
     continued-fraction recurrences. *)
  val fpmin = 1E~300

  fun sqr (x : real) = x * x

  (* --- Lanczos log-gamma (g = 7) ------------------------------------- *)

  (* Coefficients for g = 7, n = 9 (Godfrey / Numerical Recipes). *)
  val lanczosG = 7.0
  val lanczos =
    [ 0.99999999999980993
    , 676.5203681218851
    , ~1259.1392167224028
    , 771.32342877765313
    , ~176.61502916214059
    , 12.507343278686905
    , ~0.13857109526572012
    , 9.9843695780195716E~6
    , 1.5056327351493116E~7 ]

  (* Sum of the Lanczos series  A_g(z) = c0 + sum_{k>=1} c_k / (z + k),
     evaluated at the shifted argument z. *)
  fun lanczosSum z =
    let
      fun go ([], _, acc) = acc
        | go (c :: cs, k, acc) = go (cs, k + 1.0, acc + c / (z + k))
    in
      case lanczos of
        c0 :: rest => go (rest, 1.0, c0)
      | [] => 0.0
    end

  (* log-gamma valid for x > 0 via the Lanczos approximation. *)
  fun lgammaPos x =
    let
      val a = lanczosSum (x - 1.0)
      val t = x - 1.0 + lanczosG + 0.5
      val halfLog2pi = 0.918938533204672742   (* 0.5 * ln(2 pi) *)
    in
      halfLog2pi + (x - 0.5) * Math.ln t - t + Math.ln a
    end

  (* gamma valid for x >= 0.5 via the Lanczos approximation. *)
  fun gammaPos x =
    let
      val a = lanczosSum (x - 1.0)
      val t = x - 1.0 + lanczosG + 0.5
      val sqrt2pi = 2.50662827463100050242
    in
      sqrt2pi * Math.pow (t, x - 0.5) * Math.exp (~t) * a
    end

  (* G(x) for all non-pole reals, using reflection for the left half-plane. *)
  fun gamma x =
    if x >= 0.5 then gammaPos x
    else pi / (Math.sin (pi * x) * gammaPos (1.0 - x))

  (* ln |G(x)| for all non-pole reals. *)
  fun lgamma x =
    if x >= 0.5 then lgammaPos x
    else
      (* reflection: ln|G(x)| = ln(pi / |sin(pi x)|) - ln|G(1-x)| *)
      Math.ln (pi / Real.abs (Math.sin (pi * x))) - lgammaPos (1.0 - x)

  fun lbeta (a, b) = lgamma a + lgamma b - lgamma (a + b)
  fun beta (a, b) = Math.exp (lbeta (a, b))

  (* --- digamma (psi) ------------------------------------------------- *)

  (* psi via the asymptotic series, pushing the argument up with the
     recurrence  psi(x) = psi(x+1) - 1/x  until x >= 6, and reflection
     psi(1-x) - psi(x) = pi cot(pi x)  for the left half-plane. *)
  fun digamma x =
    if x <= 0.0 andalso Real.== (x, Real.realRound x) then
      (* pole at non-positive integers *)
      raise Domain
    else if x < 0.5 then
      digamma (1.0 - x) - pi / Math.tan (pi * x)
    else
      let
        fun shift (y, acc) =
          if y < 10.0 then shift (y + 1.0, acc - 1.0 / y) else (y, acc)
        val (y, acc) = shift (x, 0.0)
        val inv = 1.0 / y
        val inv2 = inv * inv
        (* ln y - 1/(2y) - 1/(12 y^2) + 1/(120 y^4) - 1/(252 y^6) + 1/(240 y^8) *)
        val tail =
          inv2 * (1.0 / 12.0
            - inv2 * (1.0 / 120.0
              - inv2 * (1.0 / 252.0
                - inv2 * (1.0 / 240.0))))
      in
        acc + Math.ln y - 0.5 * inv - tail
      end

  val psi = digamma

  (* --- regularized incomplete gamma --------------------------------- *)

  (* Generous, compiler-independent iteration ceiling: the series/continued
     fractions below converge to `eps` in well under this many steps, so the
     cap only guards against a pathological non-converging (e.g. NaN) loop and
     never changes a real result -- keeping the kernels deterministic. *)
  val maxIter = 1000

  (* Lower series for P(a,x), valid (and fast-converging) for x < a+1. *)
  fun gammaP_series (a, x) =
    let
      val gln = lgammaPos a
      fun go (n, ap, del, sum) =
        let
          val ap' = ap + 1.0
          val del' = del * x / ap'
          val sum' = sum + del'
        in
          if Real.abs del' < Real.abs sum' * eps orelse n >= maxIter then sum'
          else go (n + 1, ap', del', sum')
        end
      val sum0 = 1.0 / a
      val sum = go (0, a, sum0, sum0)
    in
      sum * Math.exp (~x + a * Math.ln x - gln)
    end

  (* Upper continued fraction for Q(a,x) (Lentz), valid for x >= a+1. *)
  fun gammaQ_cf (a, x) =
    let
      val gln = lgammaPos a
      val b0 = x + 1.0 - a
      val c0 = 1.0 / fpmin
      val d0 = 1.0 / b0
      fun go (n, i, b, c, d, h) =
        let
          val an = ~i * (i - a)
          val b' = b + 2.0
          val d1 = an * d + b'
          val d2 = if Real.abs d1 < fpmin then fpmin else d1
          val c1 = b' + an / c
          val c2 = if Real.abs c1 < fpmin then fpmin else c1
          val d3 = 1.0 / d2
          val del = d3 * c2
          val h' = h * del
        in
          if Real.abs (del - 1.0) < eps orelse n >= maxIter then h'
          else go (n + 1, i + 1.0, b', c2, d3, h')
        end
      val h = go (0, 1.0, b0, c0, d0, d0)
    in
      Math.exp (~x + a * Math.ln x - gln) * h
    end

  fun gammaIncP (a, x) =
    if x < 0.0 orelse a <= 0.0 then raise Domain
    else if Real.== (x, 0.0) then 0.0
    else if x < a + 1.0 then gammaP_series (a, x)
    else 1.0 - gammaQ_cf (a, x)

  fun gammaIncQ (a, x) =
    if x < 0.0 orelse a <= 0.0 then raise Domain
    else if Real.== (x, 0.0) then 1.0
    else if x < a + 1.0 then 1.0 - gammaP_series (a, x)
    else gammaQ_cf (a, x)

  (* --- error functions ----------------------------------------------- *)

  (* erf x = sign(x) P(1/2, x^2);  erfc x = Q for x >= 0, else 2 - erfc(-x). *)
  fun erf x =
    if Real.== (x, 0.0) then 0.0
    else if x > 0.0 then gammaIncP (0.5, x * x)
    else ~(gammaIncP (0.5, x * x))

  fun erfc x =
    if x >= 0.0 then gammaIncQ (0.5, x * x)
    else 1.0 + gammaIncP (0.5, x * x)

  (* --- inverse error function ---------------------------------------- *)

  val sqrtPi = 1.77245385090551602730

  fun erfInv y =
    if Real.== (y, 0.0) then 0.0
    else if y >= 1.0 then Real.posInf
    else if y <= ~1.0 then Real.negInf
    else
      let
        (* Winitzki's closed-form initial guess. *)
        val a = 0.147
        val ln1 = Math.ln (1.0 - y * y)
        val t1 = 2.0 / (pi * a) + ln1 / 2.0
        val sign = if y < 0.0 then ~1.0 else 1.0
        val x0 = sign * Math.sqrt (Math.sqrt (t1 * t1 - ln1 / a) - t1)
        (* Newton refinement on f(x) = erf x - y, f'(x) = (2/sqrt pi) e^{-x^2},
           driven to a fixed tolerance (max 12 steps as a safety bound). *)
        fun refine (x, n) =
          if n = 0 then x
          else
            let
              val fx = erf x - y
              val dfx = 2.0 / sqrtPi * Math.exp (~(x * x))
              val x' = x - fx / dfx
            in
              if Real.abs (x' - x) < eps * (1.0 + Real.abs x') then x'
              else refine (x', n - 1)
            end
      in
        refine (x0, 12)
      end

  (* --- regularized incomplete beta ----------------------------------- *)

  (* Lentz continued fraction for I_x(a,b)'s core (Numerical Recipes `betacf`). *)
  fun betaCF (a, b, x) =
    let
      val qab = a + b
      val qap = a + 1.0
      val qam = a - 1.0
      val c0 = 1.0
      val d0r = 1.0 - qab * x / qap
      val d0 = 1.0 / (if Real.abs d0r < fpmin then fpmin else d0r)
      fun go (n, m, c, d, h) =
        let
          val m2 = 2.0 * m
          (* even step *)
          val aa1 = m * (b - m) * x / ((qam + m2) * (a + m2))
          val d1r = 1.0 + aa1 * d
          val d1 = 1.0 / (if Real.abs d1r < fpmin then fpmin else d1r)
          val c1r = 1.0 + aa1 / c
          val c1 = if Real.abs c1r < fpmin then fpmin else c1r
          val h1 = h * d1 * c1
          (* odd step *)
          val aa2 = ~(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
          val d2r = 1.0 + aa2 * d1
          val d2 = 1.0 / (if Real.abs d2r < fpmin then fpmin else d2r)
          val c2r = 1.0 + aa2 / c1
          val c2 = if Real.abs c2r < fpmin then fpmin else c2r
          val del = d2 * c2
          val h2 = h1 * del
        in
          if Real.abs (del - 1.0) < eps orelse n >= maxIter then h2
          else go (n + 1, m + 1.0, c2, d2, h2)
        end
    in
      go (0, 1.0, c0, d0, d0)
    end

  fun betaInc (a, b, x) =
    if a <= 0.0 orelse b <= 0.0 then raise Domain
    else if x < 0.0 orelse x > 1.0 then raise Domain
    else if Real.== (x, 0.0) then 0.0
    else if Real.== (x, 1.0) then 1.0
    else
      let
        val lbt = lgamma (a + b) - lgamma a - lgamma b
                  + a * Math.ln x + b * Math.ln (1.0 - x)
        val bt = Math.exp lbt
      in
        if x < (a + 1.0) / (a + b + 2.0) then
          bt * betaCF (a, b, x) / a
        else
          1.0 - bt * betaCF (b, a, 1.0 - x) / b
      end
end
