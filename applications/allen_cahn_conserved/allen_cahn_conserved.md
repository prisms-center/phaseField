# PRISMS PhaseField: Globally Conserved Allen-Cahn Dynamics

This application performs a phase field simulation of Allen-Cahn dynamics subject to global (as opposed to local) conservation. Global conservation implies that $\int_\Omega c \,\mathrm{d}V$ is constant in time.

Consider a free energy expression of the form:

$$
\begin{equation}
F[c] = \int_{\Omega} f(c) + \frac{\kappa}{2} \nabla c \cdot \nabla c \,\mathrm{d} V,
\end{equation}
$$

where $f(c)$ is the chemical free energy density, $\kappa$ is the gradient coefficient, and $c$ is the structural order parameter, which subject to a global constraint

$$
\begin{equation}
\int_{\Omega} c \,\mathrm{d} V = \mathrm{constant}.
\end{equation}
$$

The kinetics is given by the Allen-Cahn equation, which is the gradient flow of the functional $F[c]$:

$$
\begin{equation}
\frac{\partial c}{\partial t} = - M \mu,
\end{equation}
$$

where $M$ is the mobility and $\mu$ is the chemical potential defined as

$$
\begin{equation}
\mu = \frac{\delta F}{\delta c} = f_{,c}(c) - \kappa \nabla^2 c.
\end{equation}
$$

To conserve $c$, we can introduce a global Lagrange multiplier $\lambda(t)$

$$
\begin{equation}
\frac{\partial c}{\partial t} = - M (\mu - \lambda(t)).
\end{equation}
$$

The Lagrange multiplier can be solved from the conservation of $c$ (Eq. 2) as

$$
\begin{equation}
\lambda(t) = \frac{1}{V} \int_{\Omega} \mu \,\mathrm{d} V, \quad V = \int_{\Omega} \,\mathrm{d} V.
\end{equation}
$$

The governing equations are given by Eqs. 5, 4, and 6.


## Time discretization
Considering forward Euler explicit time stepping, we have the time discretized kinetics equations:

$$
\begin{equation}
c^{n+1} = c^{n} - \Delta t M(\mu^n-\lambda^n)
\end{equation}
$$

and

$$
\begin{equation}
\mu^{n} = f_{,c}^n - \kappa \nabla^2 c^n.
\end{equation}
$$

## Weak formulation
In the weak formulation, considering an arbitrary variation $w$, the above equation can be expressed as a residual equations:

$$
\begin{equation}
\int_{\Omega} w c^{n+1} \,\mathrm{d}V = \int_{\Omega} w \left( c^{n} - \Delta t M(\mu^n-\lambda^n) \right) \,\mathrm{d}V
\end{equation}
$$

$$
\begin{equation}
r_{c} = c^{n} - \Delta t M(\mu^n-\lambda^n)
\end{equation}
$$

and

$$
\begin{align}
\int_{\Omega} w \mu^{n} \,\mathrm{d}V &= \int_{\Omega} w f_{,c}^n - w \kappa \nabla^2 c^{n} \,\mathrm{d}V \\
&= \int_{\Omega} w ~ (f_{,c}^{n}) + \nabla w \cdot (\kappa \nabla c^{n}) \,\mathrm{d}V,
\end{align}
$$

$$
\begin{equation}
r_{\mu} = f_{,\eta}^{n}
\end{equation}
$$

$$
\begin{equation}
r_{\mu x} = \kappa \nabla \eta^{n}
\end{equation}
$$

where the Lagrange multiplier $\lambda^n$ for time step $n$ is calculated as

$$
\begin{align}
\lambda^n& = \frac{1}{V}\int_\Omega \mu^n \,\mathrm{d}V.
\end{align}
$$

The above values of $r_{c}$, $r_{\mu}$, and $r_{\mu x}$ are used in `compute_rhs()` in the following source file:
`applications/allen_cahn_conserved/custom_pde.cc`
