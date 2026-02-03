= Optimal control strategy <supp-mat-optcontrol>
The optimal control problem is stated as follows. Find the optimal control strategy $alpha(t)$ that minimises the following cost function $J$:
$
  J = integral_(t_0)^(t_f) dot(nu) dif t
$
subject to the system dynamics:
$
     dot(nu) & = -C_D nu^2 - sin gamma \
  dot(gamma) & = C_L nu - (cos gamma) / nu
$
the final time, $t_f$, is left as a free variable.

== Derivation of optimal control strategy
The Hamiltonian, $cal(H)$, for this problem is formed by introducing the costate variables, $p_1, p_2$:
$
  cal(H) & = dot(nu) + p_1 dot(nu) + p_2 dot(gamma) \
         & = (1 + p_1)(-C_D nu^2 -sin gamma) + p_2 (C_L nu - cos gamma \/ nu)
$
Using Pontryagin's maximum principle, the necessary conditions for the optimal solution are @kirkOptimalControlTheory2004:
$
    (partial cal(H)) / (partial p_1) & = dot(nu) \
    (partial cal(H)) / (partial p_2) & = dot(gamma) \
     (partial cal(H)) / (partial nu) & = -dot(p_1) \
  (partial cal(H)) / (partial gamma) & = -dot(p_2) \
  (partial cal(H)) / (partial alpha) & = 0 \
$
along with the boundary conditions:
$
     p_1(t_f) & = 0 \
     p_2(t_f) & = 0 \
  cal(H)(t_f) & = 0
$
which result from the final time being free. This results in a system of ordinary differential equations given by:
$
     dot(nu) & = -C_D nu^2 - sin gamma \
  dot(gamma) & = C_L nu - cos gamma \/ nu \
    dot(p_1) & = -(1+p_1)(-2 C_D nu) - p_2 (C_L + cos gamma \/ nu^2) \
    dot(p_2) & = (1+p_1) cos gamma - p_2 sin gamma \/ nu \
$
with the constants of integration:
$
     nu(t_0) & = nu_0 \
  gamma(t_0) & = gamma_0 \
    p_1(t_f) & = 0 \
    p_2(t_f) & = 0
$
The optimal control law is given by:
$
  (partial cal(H)) / (partial alpha) = -(1+p_1) nu^2 (partial C_D) / (partial alpha) + p_2 nu (partial C_L) / (partial alpha) = 0
$ <eq-optcontrol-costate>
which can be solved explicitly for $alpha$ in terms of the state and costate variables for the aerodynamics model used here (@eq-clcd):
$
  cot 2alpha = (1+p_1)/p_2 C_(D"max")/(2C_(L"max")) nu
$
Conventional techniques for solving boundary value problems can then be used to determine the optimal solution. However, this problem can more easily solved as an initial value problem by choosing a touchdown speed. For a given touchdown speed $nu_("TD")$, the resulting value of $gamma$ can be calculated such that it gives $dot(nu) = 0$ (i.e. the minimum speed has been reached). To do so however requires knowing the value of $C_D$ (and therefore $alpha$) at touchdown. This can be found by considering @eq-optcontrol-costate evaluated at $t_f$:
$
  (partial cal(H)) / (partial alpha) (t_f) = -nu^2 (partial C_D) / (partial alpha) = 0
$
as $p_1(t_f) = p_2(t_f) = 0$. Thus choosing the value of $alpha$ which maximises $C_D$ is the optimal at $t_f$; optimal trajectories touchdown at maximum $C_D$. As $nu$, $gamma$, $p_1$, and $p_2$ are all known at $t_f$ the solution can be integrated backwards in time to an arbitrary point to generate the solution.

Alternatively, the optimal control law can be derived without the inclusion of the costate variables, as given in the text (@eq-opt-control-alpha). If the final time is left free, as is the case here, then $cal(H) = 0$ for $t in [t_0, t_f]$ along the optimal solution @kirkOptimalControlTheory2004, therefore:
$
  (1 + p_1) dot(nu) + p_2 dot(gamma) = 0
$
or equivalently:
$
  dot(gamma) / dot(nu) = - (1 + p_1) / p_2
$<eq-opt-ham>
Rewriting @eq-optcontrol-costate as:
$
  (partial C_L \/ partial alpha) / (partial C_D \/ partial alpha) = (1 + p_1) /p_2 nu
$
and combining with @eq-opt-ham:
$
  dot(gamma) / dot(nu) = -1 / nu (partial C_L \/ partial alpha) / (partial C_D \/ partial alpha)
$
or:
$
  (partial C_D) / (partial alpha) nu dot(gamma) + (partial C_L) / (partial alpha) dot(nu) = 0
$ <eq-opt-control-pmp>
which is the optimal control law as given in @eq-opt-control-alpha. This equation can then be solved for $alpha$ as a function of just the state variables $nu$ and $gamma$. Here, this equation is solved for $alpha$ numerically using the ITP method, which has similar robustness to the bisection method but with improved convergence. The ITP method requires an interval with opposite signs of the function at each end, ensuring the solution is enclosed in the interval. This interval is given by $alpha_("min")$, which is the minimum value of $alpha$ for which $dot(nu) < 0$ and $pi\/2$.

Alternatively the optimal control law can be derived from assuming that the optimal strategy is to maximise the increase in flight path angle per decrease in speed, or $dot(gamma) / dot(nu)$, as given in the main text. This ratio can be expressed as:
$
  dot(gamma) / dot(nu) = (C_L nu - cos gamma \/ nu) / (-C_D nu^2 - sin gamma)
$
taking the derivative with respect to $alpha$:
$
  (partial) / (partial alpha) dot(gamma) / dot(nu) = (dot(gamma) (partial dot(nu)) / (partial alpha) - dot(nu) (partial dot(gamma)) / (partial alpha)) / dot(nu)^2
$
the value of $alpha$ which minimises $dot(gamma) / dot(nu)$ must occur when the above equation equals 0:
$
  dot(gamma) (partial dot(nu)) / (partial alpha) - dot(nu) (partial dot(gamma)) / (partial alpha) = 0
$
which can simplify to:
$
  dot(gamma) nu (partial C_D) / (partial alpha) + dot(nu) (partial C_L) / (partial alpha) = 0
$<eq-opt-control2>
which is the same as derived above using the necessary conditions for optimality (@eq-opt-control-pmp).

== Solutions to optimal control strategy
The optimal control problem can thus be solved in three distinct ways

- solve the augmented system of state and co-state variables with the optimal control given as the solution to @eq-optcontrol-costate. The resulting system can either be solved as a two-point boundary value problem (BVP), with the state known and the initial time and the co-state known at the final time, or, by specifying the touchdown speed, an initial value problem (IVP) solved backwards in time from the known state and co-state at touchdown.

- solve the original system (no co-state) with the optimal control strategy given by numerically solving @eq-opt-control-pmp/@eq-opt-control2 at each time step.

- solve the optimal control problem directly as a non-linear programming problem using direct collocation as described in @sec-methods-optcontrol.

The three different methods produce identical solutions within the error of the numerical methods (@fig-opt-control-compare).


#figure(
  image("../../figs/optcontrol_comparison.pdf"),
  caption: [Comparison between the different methods to solve the optimal control problem. BVP: solved as augmented system of state and co-state variables with split boundary conditions. Here BVP is solved as an initial value problem backwards from touchdown. IVP: solves original system numerically solving for $alpha$ at each time step. DC: direct solution to optimal control problem using direct collocation. Trajectories solved for a touchdown speed of 0.1.],
) <fig-opt-control-compare>
