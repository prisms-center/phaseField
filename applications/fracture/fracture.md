# PRISMS PhaseField: Phase-Field Fracture

This application implements a phase-field model of brittle fracture, closely following
Hossain et al. [1], with the phase field evolved by a forward-Euler time iteration rather
than solved as a stationary point of the energy, in the spirit of Kuhn and Müller [2]. Crack
growth is driven by a **surfing boundary condition**: a mode-I asymptotic displacement field
that translates through the domain at constant velocity, holding a steady stress-intensity
factor near a virtual crack tip.

## Governing equations

Consider a total free energy of the form

$$
\begin{equation}
\Pi(u,n) = \int_{\Omega} E(\mathbf{x})\,h(n)\,\Psi(u) ~dV + \int_{\Omega} G_c(\mathbf{x}) \left( \frac{3}{8\ell} n + \frac{3\ell}{8} |\nabla n|^2 \right) dV - \int_{\partial \Omega} u \cdot t ~dS
\end{equation}
$$

where $u$ is the displacement and $n \in [0,1]$ is a scalar phase field describing damage
($n=0$: intact, $n=1$: fully broken). $\Psi(u) = \frac{1}{2}\varepsilon(u):C_0:\varepsilon(u)$
is the elastic strain energy density of the undamaged base material ($\varepsilon(u) =
\frac{1}{2}(\nabla u + \nabla u^T)$, $C_0$ the base elasticity tensor), and $h(n)=(1-n)^2$ is
the quadratic degradation function. $E(\mathbf{x}) \in [0,1]$ and $G_c(\mathbf{x})$ are
dimensionless masks (fields `Ex`, `Gx` in code) that let the base stiffness $C_0$ and baseline
toughness $G_c^0$ vary in space, for studies of heterogeneous media [1]; they are assigned
once through the initial condition and held fixed thereafter. $t = \sigma \cdot \hat{n}$ is
the surface traction; body forces are zero.

## Variational treatment

**Mechanical equilibrium.** Holding $n$ fixed and varying $u \to u + \alpha w$:

$$
\begin{align}
\delta_u \Pi &= \int_{\Omega} \nabla w : \sigma ~dV - \int_{\partial \Omega} w \cdot t ~dS \\
\sigma &= E(\mathbf{x})\,h(n)\,C_0:\varepsilon(u)
\end{align}
$$

Loading is applied through Dirichlet conditions on $u$ (see Surfing boundary condition
below), so $\delta_u \Pi = 0$ reduces to the weak form $R_u(w) = \int_{\Omega} \nabla w : \sigma ~dV = 0$,
i.e. $\nabla \cdot \sigma = 0$.

**Phase-field evolution.** $n$ evolves by a non-conserved, $L^2$ gradient-flow kinetic law
with mobility $M_n$ [2]; this time dependence is a numerical regularization rather than a
physical rate law [4], so loading rates should be kept slow enough that $M_n$ does not affect
the solution. Varying $n \to n + \alpha v$ at fixed $u$:

$$
\begin{align}
\delta_n \Pi &= \int_{\Omega} E(\mathbf{x})\,h_{,n}\,\Psi(u)\,v ~dV + \int_{\Omega} G_c(\mathbf{x}) \left( \frac{3}{8\ell} v + \frac{3\ell}{4} \nabla n \cdot \nabla v \right) dV, \qquad h_{,n} = 2(n-1)
\end{align}
$$

giving the evolution law $\int_{\Omega} \dot{n}~v ~dV = -M_n \delta_n \Pi[v]$, or in strong form

$$
\begin{equation}
\dot{n} = -M_n \left[ 2(n-1)\,E(\mathbf{x})\Psi(u) + G_c(\mathbf{x}) \frac{3}{8\ell} - G_c(\mathbf{x}) \frac{3\ell}{4} \nabla^2 n \right]
\end{equation}
$$

subject to the **irreversibility constraint** $\dot n \ge 0$ (damage cannot heal) and the box
constraint $n \le 1$, both enforced numerically by clamping the explicit update at each time
step (see Time discretization). Splitting $\Psi(u)$ into tensile/compressive parts to prevent
crack growth under pure compression [3] is not implemented; boundary conditions should be
chosen to induce predominantly tensile stress states.

## Kinetics

Mechanical equilibrium is elliptic and is re-solved every time step for the current damage
field, while the phase field evolves according to a parabolic rate law:

$$
\begin{align}
\nabla \cdot \sigma &= 0 \\
\sigma &= E(\mathbf{x})\,h(n)\,C_0:\varepsilon(u)
\end{align}
$$

$$
\begin{align}
\frac{\partial n}{\partial t} &= -M_n \left[ 2(n-1)\,E(\mathbf{x})\Psi(u) + G_c(\mathbf{x}) \frac{3}{8\ell} - G_c(\mathbf{x}) \frac{3\ell}{4} \nabla^2 n \right]
\end{align}
$$

## Time discretization

$\dot n$ is computed once per step as an auxiliary field (`dndt` in code) so that $\Psi(u)$
is evaluated only once and reused, then $n$ is advanced by forward Euler:

$$
\begin{align}
\dot{n}^{n-1} &= -M_n \left[ 2(n^{n-1}-1)\,E(\mathbf{x})\Psi(u^{n-1}) + G_c(\mathbf{x}) \frac{3}{8\ell} - G_c(\mathbf{x}) \frac{3\ell}{4} \nabla^2 n^{n-1} \right]
\end{align}
$$

$$
\begin{align}
n^{n} &= n^{n-1} + \Delta t\, \dot{n}^{n-1}
\end{align}
$$

clamped afterward so that the update satisfies $n^{n} \ge n^{n-1}$ (irreversibility) and
$n^{n} \le 1$ (box constraint).

## Weak formulation

For the auxiliary field $\dot n$, with arbitrary variation $w$:

$$
\begin{align}
\int_{\Omega} w\, \dot{n}^{n-1} ~dV &= \int_{\Omega} w\, \mathrm{RHS}_{dndt} + \nabla w \cdot \mathrm{RHS}_{dndtx} ~dV
\end{align}
$$

$$
\begin{align}
\mathrm{RHS}_{dndt} &= -M_n \left[ 2(n^{n-1}-1)\,E(\mathbf{x})\Psi(u^{n-1}) + G_c(\mathbf{x}) \frac{3}{8\ell} \right]
\end{align}
$$

$$
\begin{align}
\mathrm{RHS}_{dndtx} &= -M_n\, G_c(\mathbf{x}) \frac{3\ell}{4} \nabla n^{n-1}
\end{align}
$$

For the phase field $n$ itself, no gradient term is needed (the Laplacian has already been
folded into $\dot n$ above) so the update is a pointwise sum:

$$
\begin{align}
\int_{\Omega} w\, n^{n} ~dV &= \int_{\Omega} w\, \mathrm{RHS}_{n} ~dV, \qquad \mathrm{RHS}_{n} = n^{n-1} + \Delta t\, \dot{n}^{n-1}
\end{align}
$$

For the displacement field $u$, loading enters only through the Dirichlet boundary condition,
so the load-vector side of the linear solve carries no source term ($\mathrm{RHS}_{u} = 0$,
$\mathrm{RHS}_{ux} = 0$); the matrix-free CG solve instead assembles, at each iteration, the action of
the degraded tangent stiffness on a trial update $\Delta u$:

$$
\begin{align}
\mathrm{LHS}_{ux} &= E(\mathbf{x})\,h(n)\,C_0:\varepsilon(\Delta u)
\end{align}
$$

The above expressions of $\mathrm{RHS}_{dndt}$, $\mathrm{RHS}_{dndtx}$, $\mathrm{RHS}_{n}$, and $\mathrm{LHS}_{ux}$ define the
code written in:
`custom_pde.h`

## Surfing boundary condition

Crack growth is driven entirely by a Dirichlet condition on $u$: the leading-order,
plane-strain mode-I asymptotic displacement field for a semi-infinite crack in an infinite
body [5], evaluated in a frame centered on a virtual crack tip that translates at constant
velocity $v_{\text{nom}}$: the "surfing" boundary condition of Hossain et al. [1]:

$$
\begin{align}
u_x &= \frac{K_I^{\text{nom}}}{2\mu}\sqrt{\frac{r}{2\pi}} (\kappa - \cos\theta)\cos\frac{\theta}{2} \\
u_y &= \frac{K_I^{\text{nom}}}{2\mu}\sqrt{\frac{r}{2\pi}} (\kappa - \cos\theta)\sin\frac{\theta}{2}
\end{align}
$$

where $r,\theta$ are polar coordinates centered on the moving tip
$(v_{\text{nom}}t + c_\ell,\, y_{\text{tip}})$, $\mu,\lambda$ are the Lamé parameters recovered
from $C_0$, $\nu = \lambda/[2(\lambda+\mu)]$, and $\kappa = 3-4\nu$ (plane strain). The crack
should propagate through some combination of $K_I^{\text{nom}}$ exceeding the critical
stress intensity factor $K_{Ic} = \sqrt{G_c E'}$ (with $E'=E$ for plane stress, $E'=E/(1-\nu^2)$
for plane strain) and the actual crack tip lagging the imposed tip location; $v_{\text{nom}}$
should be slow enough that the phase field stays close to equilibrium.

## Initial conditions

$n$ is seeded at $t=0$ with a straight crack of length $c_\ell$ along the domain mid-height,
using the analytical AT1 energy-minimizing profile

$$
\begin{equation}
n_0(\mathbf{x}) = \left[ \left(1 - \frac{d(\mathbf{x})}{2\ell}\right)_+ \right]^2
\end{equation}
$$

where $d(\mathbf{x})$ is the distance from $\mathbf{x}$ to the seed crack, so the crack starts
at its steady-state regularized width and does not need to re-equilibrate before propagation
begins.

## References

[1] M.Z. Hossain, C.-J. Hsueh, B. Bourdin, and K. Bhattacharya, Effective toughness of
heterogeneous media. *J. Mech. Phys. Solids* **71** (2014) 15–32. doi: 10.1016/j.jmps.2014.06.002.

[2] C. Kuhn and R. Müller, A continuum phase field model for fracture. *Eng. Fract. Mech.* **77**
(2010) 3625–3634. doi: 10.1016/j.engfracmech.2010.08.009.

[3] L. De Lorenzis and C. Maurini, Nucleation under multi-axial loading in variational phase-field
models of brittle fracture. *Int. J. Fract.* doi: 10.1007/s10704-021-00555-6.

[4] C. Miehe, F. Welschinger, and M. Hofacker, Thermodynamically consistent phase-field models of
fracture: Variational principles and multi-field FE implementations. *Int. J. Numer. Methods Eng.*
**83** (2010) 1273–1311. doi: 10.1002/nme.2861.

[5] A. T. Zehnder, *Fracture Mechanics*. Springer Dordrecht (2012). doi: 10.1007/978-94-007-2595-9.
