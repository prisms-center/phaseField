## PRISMS-PF: Advection Equation
The advection equation is given by:

$$
\begin{align}
\frac{\partial \eta}{\partial t} = -\textbf{u}\cdot\nabla\eta
\end{align}
$$

where $\eta$ is the scalar conserved quantity and $\textbf{u}$ is a constant vector field.

## Time discretization

Considering backward Euler time stepping, we have the time discretized kinetics equation:

$$
\begin{align}
\eta^{n+1} = \eta^{n} - \Delta t (\textbf{u}\cdot\nabla\eta^{n+1})
\end{align}
$$

PRISMS-PF makes use of the Newton's method for linear and nonlinear iterations, so we can substitute $\eta^{n+1}=\Delta \eta^{n+1}+n_0^{n+1}$.

$$
\begin{align}
\Delta \eta^{n+1} + \Delta t (\textbf{u}\cdot\nabla(\Delta \eta^{n+1}))= \eta^{n} - \eta_0^{n+1} - \Delta t (\textbf{u}\cdot\nabla\eta_0^{n+1})
\end{align}
$$

## Weak formulation

In the weak formulation, considering an arbitrary variation $w$, the above equation can be expressed as a residual equation:

$$
\begin{align}
\int_{\Omega} w [\Delta \eta^{n+1} + \Delta t (\textbf{u}\cdot\nabla(\Delta \eta^{n+1}))] ~dV &= \int_{\Omega} w [\eta^{n} - \eta_0^{n+1} - \Delta t (\textbf{u}\cdot\nabla\eta_0^{n+1})] ~dV
\end{align}
$$