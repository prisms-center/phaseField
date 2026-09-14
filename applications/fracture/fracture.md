# Phase-Field Fracture

This example application implements a 2D phase-field model of brittle fracture, driven by a
**surfing boundary condition**: a mode-I asymptotic displacement field that translates through the
domain at a prescribed velocity, imposing a steady stress-intensity factor near a (virtual) moving
crack tip and driving crack propagation.

## Governing equations

Consider a total free energy of the form

$$
\Pi(u,n) = \int_{\Omega} g(n) \psi_e\big(\varepsilon(u)\big) dV
 + \int_{\Omega} G_c \gamma(n,\nabla n) dV
 - \int_{\partial \Omega} u \cdot t dS
$$

where $u$ is the displacement, $\varepsilon(u) = \frac{1}{2}(\nabla u + \nabla u^T)$ is the
infinitesimal strain tensor, and $n \in [0,1]$ is a scalar phase field describing the state of
damage ($n=0$: intact material, $n=1$: fully broken). The elastic strain energy density of the
undamaged material is

$$
\psi_e(\varepsilon) = \frac{1}{2} \varepsilon : C_0 : \varepsilon ,
$$

with $C_0$ the (isotropic) elasticity tensor of the base material, and $g(n) = (1-n)^2$ is a
quadratic degradation function that reduces stiffness as the material accumulates damage. The
crack surface density function

$$
\gamma(n,\nabla n) = \frac{3}{8\ell} n + \frac{3\ell}{8} |\nabla n|^2
$$

is the (AT1, linear-softening) regularized approximation to the sharp-crack surface measure, so
that $G_c\gamma$ integrates to the Griffith fracture energy $G_c$ per unit crack length as the
regularization length $\ell \to 0$. $t = \sigma \cdot \hat{n}$ is the surface traction; body forces
are assumed to be zero.

### Mechanical equilibrium

Holding $n$ fixed and considering variations of the displacement $u \to u + \alpha w$,

$$
\begin{aligned}
\delta_u \Pi &= \frac{d}{d\alpha}\Pi(u+\alpha w, n)\Big|_{\alpha=0} \\
&= \int_{\Omega} g(n) \varepsilon(w) : C_0 : \varepsilon(u) dV - \int_{\partial\Omega} w\cdot t dS \\
&= \int_{\Omega} \nabla w : \sigma dV - \int_{\partial\Omega} w\cdot t dS ,
\end{aligned}
$$

where

$$
\sigma = g(n) C_0 : \varepsilon(u)
$$

is the degraded Cauchy stress. In this application all loading is applied through a prescribed
(time-dependent) Dirichlet condition on $u$ rather than through surface tractions, so
$\delta_u \Pi = 0$ for admissible variations $w$ reduces to the weak form

$$
R_u(w) = \int_{\Omega} \nabla w : \sigma dV = 0 \qquad \Longleftrightarrow \qquad \nabla\cdot\sigma = 0 \quad \text{in} \quad \Omega ,
$$

the usual quasi-static equilibrium equation for a linear elastic solid with spatially and
temporally varying (damage-degraded) stiffness.

### Phase-field evolution

Unlike $u$, the damage field $n$ is not required to satisfy stationarity of $\Pi$ at every
instant; instead it evolves by a non-conserved, $L^2$ gradient-flow (Ginzburg–Landau-type) kinetic
law with mobility $M_n$, in the spirit of the phase-field fracture models of Kuhn and Müller [2]
and the variational formulations reviewed in [3, 4]. Taking the variation of $\Pi$ with respect to
$n \to n + \alpha v$ at fixed $u$,

$$
\begin{aligned}
\delta_n \Pi &= \frac{d}{d\alpha}\Pi(u, n+\alpha v)\Big|_{\alpha=0} \\
&= \int_{\Omega} g'(n) \psi_e v dV + \int_{\Omega} G_c\left(\frac{3}{8\ell} v + \frac{3\ell}{4} \nabla n\cdot\nabla v\right)dV ,
\end{aligned}
$$

with $g'(n) = -2(1-n) = 2(n-1)$. The evolution law is then

$$
\int_{\Omega} \dot{n} v dV = -M_n \delta_n\Pi[v] \qquad \forall v ,
$$

or, in strong form,

$$
\dot{n} = -M_n\left[ 2(n-1) \psi_e + G_c\frac{3}{8\ell} - G_c\frac{3\ell}{4}\nabla^2 n \right] .
$$

This is subject to the **irreversibility constraint** $\dot n \ge 0$ (damage cannot heal) and the
box constraint $n \le 1$ (damage cannot exceed the fully broken state), both enforced numerically
(see below) rather than analytically.

## Surfing boundary condition

Crack growth is driven entirely by a Dirichlet condition on $u$ along $\partial\Omega$: the
classical linear-elastic (Williams) mode-I asymptotic displacement field, evaluated in a frame
centered on a virtual crack tip that translates at a constant prescribed velocity $v_{\text{nom}}$.
This is the "surfing" boundary condition of Hossain, Hsueh, Bourdin, and Bhattacharya [1], designed
to hold a nominal stress-intensity factor $K_I^{\text{nom}}$ approximately constant near the
propagating crack, so that a steady-state propagation regime can be reached without needing to
resolve the far-field problem.

Let $(x_{\text{tip}}(t), y_{\text{tip}}) = \big(v_{\text{nom}} t + c_{\ell}, L_y/2\big)$ denote
the (moving) reference point, with $c_\ell$ the initial seed-crack length and $L_y$ the domain
height. Define local polar coordinates centered on this point,

$$
x = X - x_{\text{tip}}(t), \qquad y = Y - y_{\text{tip}}, \qquad r = \sqrt{x^2+y^2}, \qquad \theta = \text{atan2}(y,x) ,
$$

for a point $(X,Y) \in \partial\Omega$. The imposed displacement is the plane-strain mode-I
near-tip field (see e.g. Zehnder [5]),

$$
\begin{aligned}
u_x(x,y,t) &= \frac{K_I^{\text{nom}}}{2\mu}\sqrt{\frac{r}{2\pi}} (\kappa - \cos\theta) \cos\frac{\theta}{2} , \\
u_y(x,y,t) &= \frac{K_I^{\text{nom}}}{2\mu}\sqrt{\frac{r}{2\pi}} (\kappa - \cos\theta) \sin\frac{\theta}{2} ,
\end{aligned}
$$

where $\mu$ and $\lambda$ are the Lamé parameters recovered from the base elasticity tensor $C_0$,
$\nu = \lambda / [2(\lambda+\mu)]$ is Poisson's ratio, and $\kappa = 3-4\nu$ is the plane-strain
Kolosov constant. As $t$ increases, this field slides through the domain at velocity
$v_{\text{nom}}$, so the boundary loading looks, in the crack-tip frame, like a stationary K-field
"surfing" past the material, hence the name.

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
