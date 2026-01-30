#let ctmax = $C_(T"max")$

== Comparison with experimental data <sec-results-perching>
Observations of Harris hawks (#emph[Parabuteo unicinctus]) flying between two perches show a learned behaviour where the birds produced deep, swooping trajectories @kleinheerenbrinkOptimizationAvianPerching2022. The previous hypothesis for this behaviour was that it minimises the distance flown post-stall @kleinheerenbrinkOptimizationAvianPerching2022 @wuestAgilePerchingManeuvers2024. In this section we use our model to find optimal solutions of this perching behaviour under different cost functions and propose an alternative explanation for these swooping flight trajectories.

The present model can simulate this flight behaviour by including thrust. Thrust was added following previous work @kleinheerenbrinkOptimizationAvianPerching2022 @wuestAgilePerchingManeuvers2024 to model birds flapping their wings during the dive phase. The thrust acts parallel to the velocity modifying the dynamics (@eq-eom-nondim) to:
$
     dot(nu) & = C_T - C_D nu^2 - sin gamma \
  dot(gamma) & = C_L nu - (cos gamma) / nu
$
where $C_T$ represents the thrust normalised by bodyweight. This models the wingbeat period time averaged thrust, which is needed to understand the overall trajectory, but neglects temporal variations within each wingbeat. The thrust coefficient forms a second control input to the model and its maximum value is constrained; we also explore the effect of varying the maximum thrust coefficient. The flight is simulated beginning just after toe-off from the first perch to touchdown at the second. Trajectories are generated using direct collocation, with the initial velocity equal to that measured in the experimental data @kleinheerenbrinkOptimizationAvianPerching2022 and the touchdown position to be within 5% of the perch. The experimental data were converted into non-dimensional form to match our model using the provided morphological data @kleinheerenbrinkOptimizationAvianPerching2022. The same aerodynamics model and parameters were used as in previous sections.

The distance flown post-stall can be expressed as:
$
  integral_0^T v(t) dot f(alpha) dif t
$
where $T$ is the touchdown time and $f$ is an indicator function that returns 1 when $alpha > alpha_(upright("stall"))$ and 0 otherwise. This equation integrates the velocity over time but only when the angle of attack is larger than the stall angle, thus giving the distance flown while stalled. Minimising the velocity after the stall angle has been reached will therefore minimise the distance flown stalled. However as these perching trajectories finish with the angle of attack above stall (as larger values of $alpha$ are required to maintain the climb to the perch) this cost function will also result in small touchdown speeds. The optimal solution for a cost function minimising the distance flown stalled will therefore produce similar trajectories to solutions that use a cost function minimising touchdown speed (green and purple lines in @fig-perching\A). To separate stall effects from speed effects we firstly consider the cost function:
$
  J = integral_0^T f(alpha) dif t
$
which is the _time_ spent stalled, with no representation of speed. The solutions which minimise stall time lead to much larger touchdown speeds because larger speeds reach the landing site quicker, reducing the time spend stalled (red and purple lines in @fig-perching\A). As some combination of stall and touchdown speed likely influence flight behaviour, we introduce a new cost function, which accounts for these effects independently from one another:
$
  J = lambda dot v(T) + (1 - lambda) dot integral_0^T f(alpha) dif t
$

where $lambda$ is the weighting factor between the two; by varying the weighting factor the influence of stall effects and speed effects on the trajectory can be considered methodically (@fig-perching\B). The two components in the cost function were normalised to values between 0 and 1 by taking their maximum and minimum values from the two extreme cases (e.g. the minimum and maximum touchdown speeds were taken from the solutions with $lambda = 1$ and $lambda = 0$ respectively). The stall angle is taken to be 20 degrees, in line with previous work @wuestAgilePerchingManeuvers2024, but the results are qualitatively insensitive to this choice. The hyperbolic tangent function was used as a smooth approximation for the indicator function @wuestAgilePerchingManeuvers2024.

The resulting trajectories depend on both the available thrust and the weighting factor. For low values of #ctmax, $lambda$ had only a minor effect on the resulting trajectories (@fig-perching\C) as the space of possible trajectories was limited by the constraint on touchdown position. For values of #ctmax lower than 0.75 the model was unable to reach the perch for the given initial conditions. For larger values of #ctmax, the resulting trajectories diverged at the extremes of $lambda$ (@fig-perching\E,F). Solely, minimising stall time ($lambda = 0$) caused an initial climb before descending to the perch, whereas solely minimising touchdown speed ($lambda = 1$) resulted in deep, swooping trajectories. Intermediate values of $lambda$ generally resulted in the best match with the experimental data, indicating that some combination of minimising stall effects and touchdown speed drive this behaviour in birds. Unlike as previously proposed @kleinheerenbrinkOptimizationAvianPerching2022, solely minimising the effects of stall did not produce trajectories seen in perching birds, this is likely due to these solutions having very large touchdown speeds that would be harmful to the birds (@fig-perching\B). Instead the good match previous modelling work found to the data stems from the cost function used, minimising the distance flown stalled, which conflates the effects of touchdown speed and stall.

The swooping behaviour was largely present only when touchdown speed was included in the cost function (@fig-perching\C-F), with even very small values of $lambda$ causing the trajectories to swoop downwards compared with not including touchdown speed at all. Our previous analysis of the dynamics explains why this is the case. Swooping (combined with additional thrust) allows the bird to build up more speed, leading to steeper climb angles and ultimately lower touchdown speeds. Diving plus thrust allows greater speeds to be reached than diving alone due to the spatial constraints on the landing site (@fig-perching\F).


#figure(
  image("../figs/perching.pdf"),
  // box(width: 100%, height: 150pt, fill: luma(240), stroke: 1pt + black),
  caption: [(A) Optimal solutions for minimum touchdown speed ($lambda = 1$), post-stall time (red $lambda=0$), and post-stall distance (purple) with a $C_(T"max") = 0.8$. Dashed horizontal line on $alpha$ plot is the stall angle (20 degrees). (B) Variation in touchdown speed with $lambda$, for increasing values of #ctmax from 0.8 to 1.4. (C)-(F) flight trajectories the optimal solutions with different values of $lambda$ and #ctmax. The touchdown speeds of these solutions correspond to curves in (B). Dotted black markers show the experimental data @kleinheerenbrinkOptimizationAvianPerching2022
  ],
)
<fig-perching>
