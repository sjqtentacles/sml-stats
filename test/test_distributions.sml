(* test_distributions.sml -- pdf/cdf against published closed-form values,
   plus deterministic checks on the samplers.

   Reference values:
     N(0,1): pdf(0) = 1/sqrt(2*pi) = 0.39894228..., cdf(0) = 0.5,
             cdf(1) = 0.84134475, cdf(1.96) = 0.97500210, cdf(-1) = 0.15865525
     Binom(10,0.5): pdf(5) = 252/1024 = 0.24609375, pdf(0) = 1/1024
     Poisson(3): pdf(0) = e^-3 = 0.04978707, pdf(1) = 3*e^-3 *)

structure DistributionTests =
struct
  open Support
  structure S = Stats

  val std : S.Normal.param = { mu = 0.0, sigma = 1.0 }
  val invSqrt2pi = 1.0 / Math.sqrt (2.0 * Math.pi)

  (* tolerance for the A&S erf approximation in the normal cdf *)
  val erfTol = 1E~6

  fun run () =
    let
      val () = Harness.section "Normal pdf"
      val () = checkApprox "N(0,1) pdf(0) = 1/sqrt(2pi)"
                 (invSqrt2pi, S.Normal.pdf std 0.0)
      val () = checkApprox "N(0,1) pdf symmetric"
                 (S.Normal.pdf std 1.3, S.Normal.pdf std (~1.3))
      val () = checkApprox "N(2,3) pdf(2) = 1/(3 sqrt(2pi))"
                 (invSqrt2pi / 3.0, S.Normal.pdf { mu = 2.0, sigma = 3.0 } 2.0)

      val () = Harness.section "Normal cdf"
      val () = checkApproxTol erfTol "N(0,1) cdf(0) = 0.5"
                 (0.5, S.Normal.cdf std 0.0)
      val () = checkApproxTol erfTol "N(0,1) cdf(1) = 0.8413447"
                 (0.8413447460, S.Normal.cdf std 1.0)
      val () = checkApproxTol erfTol "N(0,1) cdf(1.96) = 0.9750021"
                 (0.9750021048, S.Normal.cdf std 1.96)
      val () = checkApproxTol erfTol "N(0,1) cdf(-1) = 0.1586553"
                 (0.1586552539, S.Normal.cdf std (~1.0))
      val () = checkApproxTol erfTol "N cdf symmetry: F(x)+F(-x)=1"
                 (1.0, S.Normal.cdf std 0.7 + S.Normal.cdf std (~0.7))

      val () = Harness.section "Binomial pdf/cdf"
      val b10 : S.Binomial.param = { n = 10, p = 0.5 }
      val () = checkApprox "Binom(10,.5) pdf(5) = 252/1024"
                 (252.0 / 1024.0, S.Binomial.pdf b10 5)
      val () = checkApprox "Binom(10,.5) pdf(0) = 1/1024"
                 (1.0 / 1024.0, S.Binomial.pdf b10 0)
      val () = checkApprox "Binom(10,.5) pdf(10) = 1/1024"
                 (1.0 / 1024.0, S.Binomial.pdf b10 10)
      val () = checkApprox "Binom pdf(k<0) = 0" (0.0, S.Binomial.pdf b10 (~1))
      val () = checkApprox "Binom pdf(k>n) = 0" (0.0, S.Binomial.pdf b10 11)
      val () = checkApprox "Binom(10,.5) cdf(10) = 1" (1.0, S.Binomial.cdf b10 10)
      val () = checkApprox "Binom(10,.5) cdf(4) = 0.376953125"
                 (386.0 / 1024.0, S.Binomial.cdf b10 4)
      (* sum_{k=0}^{n} pdf(k) = 1 *)
      val () =
        let
          val total =
            List.foldl (fn (k, acc) => acc + S.Binomial.pdf b10 k) 0.0
              (List.tabulate (11, fn i => i))
        in checkApprox "Binom pdf sums to 1" (1.0, total) end

      val () = Harness.section "Poisson pdf/cdf"
      val p3 : S.Poisson.param = { lambda = 3.0 }
      val em3 = Math.exp (~3.0)
      val () = checkApprox "Poisson(3) pdf(0) = e^-3" (em3, S.Poisson.pdf p3 0)
      val () = checkApprox "Poisson(3) pdf(1) = 3 e^-3"
                 (3.0 * em3, S.Poisson.pdf p3 1)
      val () = checkApprox "Poisson(3) pdf(2) = 4.5 e^-3"
                 (4.5 * em3, S.Poisson.pdf p3 2)
      val () = checkApprox "Poisson pdf(k<0) = 0" (0.0, S.Poisson.pdf p3 (~1))
      val () = checkApprox "Poisson(3) cdf(0) = e^-3" (em3, S.Poisson.cdf p3 0)
      (* tail sum approaches 1 *)
      val () = checkApproxTol 1E~6 "Poisson cdf(40) ~ 1" (1.0, S.Poisson.cdf p3 40)

      val () = Harness.section "samplers: determinism + range"
      (* same seed -> same draw *)
      val () =
        let
          val (z1, _) = S.Normal.sample std (seeded 0w1)
          val (z2, _) = S.Normal.sample std (seeded 0w1)
        in checkApprox "Normal sample reproducible" (z1, z2) end
      val () =
        let
          val (k1, _) = S.Binomial.sample b10 (seeded 0w7)
          val (k2, _) = S.Binomial.sample b10 (seeded 0w7)
        in Harness.checkInt "Binomial sample reproducible" (k1, k2) end
      (* binomial draws are within [0, n] *)
      val () =
        let
          fun draws (0, _, acc) = acc
            | draws (m, s, acc) =
                let val (k, s') = S.Binomial.sample b10 s
                in draws (m - 1, s', (k >= 0 andalso k <= 10) andalso acc) end
        in Harness.check "Binomial draws in [0,10]" (draws (200, seeded 0w3, true)) end
      (* poisson draws are non-negative *)
      val () =
        let
          fun draws (0, _, acc) = acc
            | draws (m, s, acc) =
                let val (k, s') = S.Poisson.sample p3 s
                in draws (m - 1, s', (k >= 0) andalso acc) end
        in Harness.check "Poisson draws >= 0" (draws (200, seeded 0w9, true)) end

      val () = Harness.section "samplers: law of large numbers"
      (* mean of many N(5,2) draws is close to 5 *)
      val () =
        let
          val nDraws = 20000
          fun loop (0, _, acc) = acc
            | loop (m, s, acc) =
                let val (z, s') = S.Normal.sample { mu = 5.0, sigma = 2.0 } s
                in loop (m - 1, s', acc + z) end
          val avg = loop (nDraws, seeded 0w12345, 0.0) / real nDraws
        in checkApproxTol 0.05 "Normal sample mean ~ 5" (5.0, avg) end
      (* mean of many Binom(10,0.3) draws is close to n*p = 3 *)
      val () =
        let
          val nDraws = 20000
          val par : S.Binomial.param = { n = 10, p = 0.3 }
          fun loop (0, _, acc) = acc
            | loop (m, s, acc) =
                let val (k, s') = S.Binomial.sample par s
                in loop (m - 1, s', acc + real k) end
          val avg = loop (nDraws, seeded 0w555, 0.0) / real nDraws
        in checkApproxTol 0.05 "Binomial sample mean ~ n*p = 3" (3.0, avg) end
      (* mean of many Poisson(4) draws is close to lambda = 4 *)
      val () =
        let
          val nDraws = 20000
          fun loop (0, _, acc) = acc
            | loop (m, s, acc) =
                let val (k, s') = S.Poisson.sample { lambda = 4.0 } s
                in loop (m - 1, s', acc + real k) end
          val avg = loop (nDraws, seeded 0w8888, 0.0) / real nDraws
        in checkApproxTol 0.06 "Poisson sample mean ~ lambda = 4" (4.0, avg) end
    in
      ()
    end
end
