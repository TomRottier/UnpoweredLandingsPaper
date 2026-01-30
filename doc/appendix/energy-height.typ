#let cd0 = $C_(D 0)$
= Calculation of energy height <supp-mat-energyheight>
Due to the presence of drag the energy height corresponding to a particular speed is not as easily as determined as for conservative systems. To estimate the energy height we assume a dive at minimum, $C_D = cd0$, at a constant, vertical dive, $gamma = -pi\/2$. This results in the following equation for the non-dimensional acceleration:
$
  dot(nu)(tau) = -cd0 nu^2 + 1
$
which has an analytical solution given by:
$
  nu(tau) = (tanh(sqrt(cd0) tau)) / sqrt(cd0)
$
given an initial speed of 0. During the dive energy is lost to the air at a rate given by $cd0 nu^3$, resulting in a total energy loss of:
$
  -integral_0^tau_nu cd0 nu^3 dif tau
$
where $tau_nu$ is the time taken to reach a speed of $nu$ from rest. The integral can be evaluated to:
$
  integral_0^tau_nu = -(tanh^2(sqrt(cd0) tau_nu) - 2log(cosh(sqrt(cd0) tau_nu))) / (2cd0)
$
The non-dimensional energy height can then be calculated from the total energy balance
$
  zeta = nu^2/2 + (tanh^2(sqrt(cd0) tau_nu) - 2log(cosh(sqrt(cd0) tau_nu))) / (2cd0)
$
this corresponds to the height above some reference level that the system would need to be at to reach a speed of $nu$ at the reference level. This assumes that the model instantaneously changes its flight path angle to begin the climb phase. In reality lift will need to be produced to alter the flight path angle, which will increase the drag and energy loss. This metric thus provides a lower bound of the actual dive height required.
