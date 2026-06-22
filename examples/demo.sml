(* demo.sml

   A tour of `sml-stats` that prints a deterministic statistical summary and a
   text histogram. The PNG histogram/regression chart is deferred to `sml-plot`
   (not built yet), so the visual asset here is an ASCII histogram plus a
   summary report -- fully reproducible across MLton and Poly/ML.

   Build and run with `make example`. *)

structure S = Stats

(* Real formatting that is byte-identical across compilers (fixed decimals;
   always includes a decimal point). *)
fun fmt k x = Real.fmt (StringCvt.FIX (SOME k)) x
(* Same, but with a leading "-" (not SML's "~") for negatives. *)
fun fmtD k x =
  let val s = fmt k x
  in if String.isPrefix "~" s then "-" ^ String.extract (s, 1, NONE) else s end
fun line s = print (s ^ "\n")

(* Draw 1000 N(50, 10) samples with a fixed seed -> a reproducible stream. *)
val nSamples = 1000
val par : S.Normal.param = { mu = 50.0, sigma = 10.0 }

fun draw (0, _, acc) = List.rev acc
  | draw (m, s, acc) =
      let val (z, s') = S.Normal.sample par s
      in draw (m - 1, s', z :: acc) end

val samples = draw (nSamples, SplitMix64.seed 0w20260621, [])

val () = line "=== sml-stats demo ============================================"
val () = line ""
val () = line ("Drew " ^ Int.toString nSamples ^ " samples from Normal(mu=50, sigma=10)")
val () = line ("  with SplitMix64 seed 0x20260621 (deterministic).")
val () = line ""
val () = line "Descriptive statistics"
val () = line ("  mean      = " ^ fmt 4 (S.mean samples))
val () = line ("  stddev    = " ^ fmt 4 (S.stddev samples))
val () = line ("  variance  = " ^ fmt 4 (S.variance samples))
val () = line ("  min       = " ^ fmt 4 (S.minimum samples))
val () = line ("  q25       = " ^ fmt 4 (S.quantile 0.25 samples))
val () = line ("  median    = " ^ fmt 4 (S.median samples))
val () = line ("  q75       = " ^ fmt 4 (S.quantile 0.75 samples))
val () = line ("  max       = " ^ fmt 4 (S.maximum samples))
val () = line ""

(* ---- ASCII histogram of the samples ---- *)
val nBins = 24
val lo = S.minimum samples
val hi = S.maximum samples
val width = (hi - lo) / real nBins

val counts = Array.array (nBins, 0)
val () =
  List.app
    (fn x =>
       let
         val idx0 = Real.floor ((x - lo) / width)
         val idx = if idx0 < 0 then 0 else if idx0 >= nBins then nBins - 1 else idx0
       in
         Array.update (counts, idx, Array.sub (counts, idx) + 1)
       end)
    samples

val maxCount = Array.foldl Int.max 0 counts
val barW = 50
fun bar c =
  let val len = if maxCount = 0 then 0 else (c * barW) div maxCount
  in CharVector.tabulate (len, fn _ => #"#") end

val () = line "Histogram (Normal samples)"
val () =
  Array.appi
    (fn (i, c) =>
       let
         val binLo = lo + real i * width
         val label = StringCvt.padLeft #" " 7 (fmt 1 binLo)
       in
         line ("  " ^ label ^ " | " ^ bar c ^ " " ^ Int.toString c)
       end)
    counts
val () = line ""

(* ---- linear regression on a noisy line y = 3x + 2 + noise ---- *)
val regPts =
  let
    fun mk (0, _, acc) = List.rev acc
      | mk (i, s, acc) =
          let
            val x = real (20 - i)
            val (noise, s') = S.Normal.sample { mu = 0.0, sigma = 2.0 } s
            val y = 3.0 * x + 2.0 + noise
          in mk (i - 1, s', (x, y) :: acc) end
  in mk (20, SplitMix64.seed 0w42, []) end

val { slope, intercept, r2 } = S.linregress regPts
val () = line "Linear regression  (true model: y = 3x + 2 + N(0,2) noise)"
val () = line ("  slope     = " ^ fmt 4 slope ^ "   (true 3.0)")
val () = line ("  intercept = " ^ fmt 4 intercept ^ "   (true 2.0)")
val () = line ("  r^2       = " ^ fmt 4 r2)
val () = line ""

(* ---- t-tests ---- *)
val groupA = [5.1, 4.9, 5.3, 5.0, 5.2, 4.8, 5.1]
val groupB = [5.6, 5.8, 5.5, 5.9, 5.7, 6.0, 5.4]
val one = S.tTestOne { data = groupA, mu0 = 5.0 }
val two = S.tTestTwo { a = groupA, b = groupB }
val () = line "One-sample t-test  (group A vs mu0 = 5.0)"
val () = line ("  t = " ^ fmt 4 (#t one) ^ "   df = " ^ fmt 1 (#df one)
               ^ "   p = " ^ fmt 4 (#pValue one))
val () = line "Two-sample t-test  (group A vs group B, pooled)"
val () = line ("  t = " ^ fmt 4 (#t two) ^ "   df = " ^ fmt 1 (#df two)
               ^ "   p = " ^ fmt 6 (#pValue two))
val () = line ""

(* ---- correlation ---- *)
val cx = [1.0, 2.0, 3.0, 4.0, 5.0]
val cy = [2.0, 4.0, 5.0, 4.0, 5.0]
val () = line "Correlation  (x = [1,2,3,4,5], y = [2,4,5,4,5])"
val () = line ("  pearson   = " ^ fmtD 4 (S.pearson (cx, cy)))
val () = line ("  spearman  = " ^ fmtD 4 (S.spearman (cx, cy)))
val () = line ""

(* ---- chi-square goodness-of-fit ---- *)
val obs = [16.0, 18.0, 16.0, 14.0, 12.0, 12.0]
val exp = List.tabulate (6, fn _ => S.sum obs / 6.0)
val chi = S.chiSquareTest (obs, exp)
val () = line "Chi-square goodness-of-fit  (observed vs uniform expected)"
val () = line ("  statistic = " ^ fmtD 4 (#statistic chi)
               ^ "   df = " ^ Int.toString (#df chi)
               ^ "   p = " ^ fmtD 4 (#pValue chi))
val () = line ""

(* ---- F-test (variance ratio) ---- *)
val ftt = S.fTest (groupA, groupB)
val () = line "F-test  (variance ratio, group A vs group B)"
val () = line ("  statistic = " ^ fmtD 4 (#statistic ftt)
               ^ "   dfn = " ^ Int.toString (#dfn ftt)
               ^ "   dfd = " ^ Int.toString (#dfd ftt)
               ^ "   p = " ^ fmtD 4 (#pValue ftt))
val () = line ""
val () = line "==============================================================="
