#let stall_speed = $nu_("stall")$
#let term_speed = $nu_("term")$
#let cdp = $C_D^prime$
#let clmax = $C_(L"max")$
#let cdmax = $C_(D"max")$
#let numin = $nu_("min")$
#let gamin = $gamma_("min")$

== Topographically constrained landings <sec-results-constraints>
The previous control strategy dealt with unconstrained flight and resulted in large horizontal and vertical displacements to land. For birds, such landings may not always be possible as spatial constraints limit the possible flight trajectories. An example is landing on the ground, where the bird cannot fly below the landing site. In addition to physical constraints, biomechanical constraints could also influence the chosen flight trajectory, such as how the velocity vector is oriented relative to the leg to enable contact with a perch @roderickBirdinspiredDynamicGrasping2021. The effects of such restrictions are explored here by implementing numerical constraints on the flight path angle.

We focus on the common situation of descending to land on the ground. Ground-landings require that $gamma <= 0$ and necessitate a different solution to the optimal control problem above (@sec-results-optcontrol) as the trajectory can no longer climb to touchdown. This increases the minimum touchdown speed that can be attained compared with unconstrained landings. From the subset of possible flight trajectories, numerical solutions for the minimum touchdown speed are found for each possible touchdown angle ranging from vertically downwards ($gamma = -pi \/ 2$) to horizontal ($gamma = 0$). All trajectories begin from a trimmed glide at $nu = 1.3$ and the flight time is fixed at 3.0. These values are arbitrary and chosen to aid comparison between trajectories, changes to them only affect the relative timings of the trajectories - the minimum speeds remain the same. As analytical expressions for the optimal control strategies cannot generally be obtained, the phase spaces cannot be broken down into the regions as done above.

The optimal trajectories gradually transitioned from a vertical touchdown ($gamma = -pi\/2$, @fig-constraints red trajectory), covering a small horizontal and a large vertical distance, to a horizontal touchdown ($gamma = 0$, @fig-constraints green trajectory), covering a large horizontal and a small vertical distance. The trajectory with the minimum touchdown speed over all touchdown angles (@fig-constraints yellow trajectory) lay in between the two extremes but more closely resembled the horizontal touchdown.

The corresponding touchdown speed for each trajectory is shown in @fig-constraints\B, along with the minimum #term_speed and #stall_speed curves (i.e. when $C_D$ and $C_L$ are maximum). The touchdown speed decreases with increasing touchdown angle up to the minimum touchdown speed #numin at a touchdown angle #gamin ($-0.59 thin "rad" approx -34 thin "degrees"$). For touchdown angles less than #gamin, trajectories are able to reach the minimum terminal speed for the given touchdown angle. For touchdown angles larger than #gamin, trajectories are no longer able to reach #term_speed and so increase up to #stall_speed for $gamma = 0$.

The horizontal landing is similar to how large, fixed-wing aircraft will land from flight. The approach maintains a roughly constant flight path angle (@fig-constraints\Cii) before performing a "flare" (increase in $alpha$) just prior to touchdown. The flare increases both lift and drag, inclining the flight path to horizontal and increasing the deceleration (@fig-constraints\Ci-ii). To minimise the touchdown speed $alpha$ increases up to the value giving #clmax (@fig-constraints\Ciii) as this corresponds to the minimum #stall_speed of the system and thus is given by $nu = 1\/sqrt(clmax)$.

The minimum touchdown speed landing is an extension of the horizontal landing. Once #stall_speed has been reached, the optimal control strategy from @sec-results-optcontrol can be used as the ground-landing constraint ($gamma <= 0$) is guaranteed when $dot(gamma) < 0$. As with unconstrained flight, minimising the touchdown speed is done by reaching #stall_speed with as large a value of $gamma$ as possible. The horizontal landing achieves this and so the minimum touchdown speed landing will follow a similar strategy as the horizontal landing reaching #stall_speed at $gamma = 0$. After this point the control follows the optimal strategy given in @sec-results-optcontrol to achieve the minimum touchdown speed possible. For shallower landings ($gamin < gamma < 0$), the trajectories follow this optimal path but terminate early at the requisite touchdown angle. Thus the same control strategy can be used for any ground-landing with a touchdown angle between 0 and #gamin. The differences in the control strategy between the horizontal and minimum touchdown speed landings up to the stall speed (@fig-constraints\Ciii) stem from the fixed touchdown time used.

For steeper landings ($-pi\/2 < gamma < gamin$), the touchdown speed increases and #stall_speed is reached with $gamma < 0$. After this point, the optimal control strategy in @sec-results-optcontrol is still followed to achieve the minimum touchdown speed. For a vertical landing, #stall_speed = 0 and #term_speed can be achieved through a constant, maximum $C_D$.

The models illustrate why unpowered ground-landings involve higher touchdown speeds than unconstrained landings, which can achieve zero speed. Birds landing on the ground must therefore either use flapping wings, or a leg, or a combination of the two to brake to zero speed. This sets a minimum requirement on the forelimb and hindlimb musculature to be capable of absorbing the residual kinetic energy. This requirement may have constrained the evolution of avian morphology, in particular the redistribution of mass from the hind- to the fore-limbs @heersWingsLegsAvian2015. Unpowered theories of flight evolution @ostromBirdFlightHow1979 @norbergEvolutionVertebrateFlight1985 assume a limited ability to flap the wings and so the legs would be required to dissipate all the remaining energy, increasing the required amount of hindlimb mass. Combining palaeontological estimates of the hindlimb energy dissipation capability with the possible touchdown speeds here would provide new insights into the locomotion that was possible under unpowered theories of flight evolution.

#figure(
  image("../figs/constraints.pdf"),
  // box(width: 100%, height: 150pt, fill: luma(240), stroke: 1pt + black),
  caption: [
    Constrained optimal landing trajectories. All trajectories are
    constrained to have $gamma lt.eq 0$ and a
    touchdown angle between 0 and $- pi \/ 2$. Trajectories shown in
    bold for a touchdown angle of 0 (red), $- pi \/ 2$ (green), and the
    touchdown angle giving minimum touchdown speed (yellow). (Ci-iii) State and
    controls for these trajectories as functions of time. The dashed horizontal lines indicate the analytical solutions
    for speed and flight path angle for the corresponding cases.
  ],
)
<fig-constraints>

// Biomechanical constraints may increase the touchdown speed even further. Other factors such as the orientation of the velocity vector relative to the leg may dictate successful landing performance @roderickBirdinspiredDynamicGrasping2021 and may shift the preferred touchdown velocity away from the minimum speed solution. This would increase the hindlimb structural requirements even further. The horizontal and minimum speed solutions both travel over three times longer horizontally than the vertical touchdown solution. In environments where there is not space for large horizontal displacements (e.g. in a dense forest environment) then the landing strategy must adapt and more closely reflect the vertical touchdown strategy, leading to an increase in touchdown speed.
