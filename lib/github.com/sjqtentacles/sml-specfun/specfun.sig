(* specfun.sig -- pure Standard ML special functions.

   A dependency-free suite of the real-valued analytic special functions that
   statistics, probability and numerical work lean on: the gamma family
   (gamma / log-gamma / beta / digamma), the error functions (erf / erfc and
   the inverse erfInv), and the regularized incomplete gamma and beta
   functions.

   Everything is built from +,-,*,/ and the Basis `Math` primitives, with
   every iterative kernel driven to a fixed convergence tolerance rather than a
   fixed wall-clock or step budget, so results are deterministic and identical
   under both MLton and Poly/ML. *)

signature SPECFUN =
sig
  (* Convergence tolerance the iterative kernels (series / continued
     fractions / Newton refinement) drive to. *)
  val eps : real

  (* --- Gamma family -------------------------------------------------- *)

  (* Euler gamma function  G(x).  Defined for all reals except the
     non-positive integers (poles), via the Lanczos approximation with
     the reflection formula for x < 1/2. *)
  val gamma : real -> real

  (* Natural log of |G(x)|.  Numerically stable for large arguments where
     `gamma` itself overflows (e.g. lgamma 100 = ln 99!). *)
  val lgamma : real -> real

  (* Beta function  B(a,b) = G(a)G(b)/G(a+b),  and its natural log. *)
  val beta  : real * real -> real
  val lbeta : real * real -> real

  (* Digamma (psi) = d/dx ln G(x) = G'(x)/G(x).  `psi` is an alias. *)
  val digamma : real -> real
  val psi     : real -> real

  (* --- Error functions ----------------------------------------------- *)

  (* erf x  = (2/sqrt pi) * integral_0^x e^{-t^2} dt,  an odd function. *)
  val erf : real -> real

  (* erfc x = 1 - erf x, computed directly (via the upper incomplete gamma)
     so it stays accurate in the far tail where 1 - erf x cancels. *)
  val erfc : real -> real

  (* Inverse error function: erfInv (erf x) = x, for y in (-1, 1). *)
  val erfInv : real -> real

  (* --- Regularized incomplete functions ------------------------------ *)

  (* Lower regularized incomplete gamma  P(a,x) = g(a,x)/G(a),  a > 0, x >= 0.
     P(a,x) + Q(a,x) = 1. *)
  val gammaIncP : real * real -> real

  (* Upper regularized incomplete gamma  Q(a,x) = G(a,x)/G(a). *)
  val gammaIncQ : real * real -> real

  (* Regularized incomplete beta  I_x(a,b),  a > 0, b > 0, 0 <= x <= 1.
     I_0 = 0, I_1 = 1. *)
  val betaInc : real * real * real -> real
end
