= Rigid body dynamics and control <supp-mat-rigidbody>
Here we validate our assumption of a point mass model being a suitable abstraction to capture the flight dynamics in birds and in particular whether allowing the angle of attack to be freely chosen is a valid assumption.

#figure(
  image("../../figs/rigid_body_bird.svg", width: 100%),
  caption: [Rigid body model.],
) <fig-rigidbody>

To do so the flight dynamics are instead modelled by a rigid body. This now makes the angle of attack determined by the body orientation and flight path angle and so can no longer be freely chosen. Instead, the control input to the model is the moment arm where the resulting aerodynamic force acts (the centre of pressure (CoP) location). In general the CoP will vary with angle of attack  and typically varies between 0.25 and 0.5 of chord. Birds however are able to translate their entire wings relative to their centre of mass through sweeping the wings forwards and backwards. Here it is assumed that movement of the CoP is dominated by wing sweep such that any effect of angle of attack can be accounted for.

The dynamics are derived identically to the point mass model in that the velocity of the centre of mass of the body is represented by its speed, $v$, and flight path angle, $gamma$. The additional dynamics relate to the body orientation, $theta$.

$
             dot(v) & = -D/m - g sin gamma \
         dot(gamma) & = L/(m v) - g / v cos gamma \
  dot.double(theta) & = r / I (L cos alpha + D sin alpha)
$
where $r$ is the location of the CoP relative to the centre of mass and $I$ is the pitch moment of inertia about the centre of mass. The aerodynamic force, and thus $r$, is assumed to lie on the body axis.

As before the model is made non-dimensional through introduction of a nominal length $hat(r) = I / (hat(v) hat(t) m) = (I rho S) / (2m^2)$ and the non-dimensional moment arm $lambda$ giving:
$
            dot(nu) & = -C_D nu^2 - sin gamma \
         dot(gamma) & = C_L nu - (cos gamma) / nu \
  dot.double(theta) & =lambda nu^2 (C_L cos alpha + C_D sin alpha)
$ <eq-rigidbody-dynamics>

The same aerodynamic model is used expressing the lift and drag coefficients as sole functions of angle of attack, where now the angle of attack is given by:
$
  alpha = theta - gamma
$ <eq-rigidbody-alpha>


== Comparison with point mass model
Comparison with the point mass model can be made by specifying the control input, $lambda$, such that the rigid body model has the same angle of attack time history as the point mass model. If the required moment arm could be realistically attained then it can be assumed that the point mass is a valid approximation.

The required moment arm can be calculated by rearranging @eq-rigidbody-dynamics for $lambda$ and substituting $dot.double(theta) = dot.double(alpha) + dot.double(gamma)$, giving:
$
  lambda = (dot.double(alpha) + dot.double(gamma)) / (nu^2(C_L cos alpha + C_D sin alpha))
$ <eq-rigidbody-required-lambda>
This equation gives the moment arm required to produce an angular acceleration of the body that results in the desired pitch angle and angle of attack.

@fig-rigidbody-comparison shows the two models behaving identically when the angle of attack is matched. For the lowest touchdown speed (0.1; blue/yellow trajectory in @fig-rigidbody-comparison), the moment arm value exceeds 500 towards the end of the trajectory, due to the low speeds and therefore low aerodynamic forces requiring extremely large moment arms to produce the necessary pitch acceleration. Such values would be unrealistic for any sized bird and therefore the angle of attack strategy could not be exactly followed. However, deviations in the angle of attack have only a minor effect on the resulting flight trajectory and touchdown speeds as the aerodynamic forces are so low at this point that they have a minimal effect on the dynamics. In @fig-rigidbody-comparison, $lambda$ has been clamped between $plus.minus 10$ and the resulting effects on touchdown speed are negligible.

Conservative estimates of a pigeon sized bird ($m = 0.4 thin "kg"$, $S = 0.04 thin "m"^2$, $I = 0.02 thin "kg"thin"m"^2$), with a $lambda = 10$ corresponds to CoP shifts on the order of $1 thin "cm"$, which is easily attainable through sweeping the wings. Further, the required moment arm shifts are so small that any additional effects not considered here would be negligible. Such effects could include: changes to the aerodynamics due to wing sweep, changes to the CoP due to angle of attack, or achieving the required translation in CoP through rotating the wings rather than directly translating them. It therefore can be concluded that a point mass model is a valid approximation to model the flight dynamics.

#figure(
  image("../../figs/rigidbody_comparison.pdf"),
  caption: [Simulations comparing the point mass model (solid lines) with the rigid body model (dashed lines) at touchdown speeds of 0.1, 0.2, and 0.3. The rigid body model is controlled by the moment arm between the centre of mass and the centre of pressure of the aerodynamic force. The moment arm is chosen such that it follows the same angle of attack time history as the point mass model. The moment arm is clamped between $plus.minus 10$. Clamping the moment arm means the desired angle of attack time history cannot be followed but has minimal effect on the resulting trajectory and touchdown speed.],
) <fig-rigidbody-comparison>
