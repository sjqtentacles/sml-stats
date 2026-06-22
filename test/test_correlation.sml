(* test_correlation.sml -- Pearson and Spearman correlation.

   Closed-form vectors (epsilon 1e-9):
     pearson([1..5],[2,4,6,8,10])  = 1   (y = 2x, perfectly linear)
     pearson([1..5],[10,8,6,4,2])  = -1
     pearson([1..5],[2,4,5,4,5]):
       mx = 3, my = 4; dx = [-2,-1,0,1,2], dy = [-2,0,1,0,1]
       Sum dxdy = 6, Sum dx^2 = 10, Sum dy^2 = 6
       => r = 6/sqrt 60 = 0.7745966692
     spearman tie example: y = [1,2,2,3,4] -> average ranks [1,2.5,2.5,4,5];
       vs x ranks [1,2,3,4,5]: Sum dxdy = 9.5, dx^2 = 10, dy^2 = 9.5
       => rho = 9.5/sqrt 95. *)

structure CorrelationTests =
struct
  open Support
  structure S = Stats

  fun run () =
    let
      val () = Harness.section "pearson correlation"
      val () = checkApprox "perfect positive r = 1"
                 (1.0, S.pearson ([1.0,2.0,3.0,4.0,5.0],
                                  [2.0,4.0,6.0,8.0,10.0]))
      val () = checkApprox "perfect negative r = -1"
                 (~1.0, S.pearson ([1.0,2.0,3.0,4.0,5.0],
                                   [10.0,8.0,6.0,4.0,2.0]))
      val () = checkApprox "r = 6/sqrt 60 ~ 0.7745966692"
                 (6.0 / Math.sqrt 60.0,
                  S.pearson ([1.0,2.0,3.0,4.0,5.0],
                             [2.0,4.0,5.0,4.0,5.0]))
      val () = checkApprox "correlation is an alias for pearson"
                 (S.pearson ([1.0,2.0,3.0,4.0,5.0],[2.0,4.0,5.0,4.0,5.0]),
                  S.correlation ([1.0,2.0,3.0,4.0,5.0],[2.0,4.0,5.0,4.0,5.0]))
      val () = Harness.checkRaises "pearson length mismatch raises"
                 (fn () => S.pearson ([1.0,2.0],[1.0]))
      val () = Harness.checkRaises "pearson constant input raises"
                 (fn () => S.pearson ([1.0,1.0,1.0],[1.0,2.0,3.0]))

      val () = Harness.section "spearman rank correlation"
      val () = checkApprox "monotone increasing rho = 1"
                 (1.0, S.spearman ([1.0,2.0,3.0,4.0,5.0],
                                   [10.0,20.0,30.0,40.0,50.0]))
      val () = checkApprox "monotone nonlinear rho = 1"
                 (1.0, S.spearman ([1.0,2.0,3.0,4.0,5.0],
                                   [1.0,4.0,9.0,16.0,25.0]))
      val () = checkApprox "reversed rho = -1"
                 (~1.0, S.spearman ([1.0,2.0,3.0,4.0,5.0],
                                    [5.0,4.0,3.0,2.0,1.0]))
      val () = checkApprox "rho with a tie via average ranks"
                 (9.5 / Math.sqrt 95.0,
                  S.spearman ([1.0,2.0,3.0,4.0,5.0],
                              [1.0,2.0,2.0,3.0,4.0]))
    in
      ()
    end
end
