## PRISMS-PF: Allen-Cahn-Reaction Dynamics

The governing partial differential equation (PDE) is:

$$
\begin{equation}
\frac{\partial c}{\partial t} = R(\eta),
\end{equation}
$$

where $c$ is concentration, $R(\eta)$ is given by the reaction kinetics, $\eta$ is the overpotential

$$
\eta =  \mu + V(t).
$$

Here $V$ is the applied voltage and $\mu$ is the chemical potential

$$
\mu = f_{,c}(c) - \kappa \nabla^2 c,
$$

where $f$ is the free energy density, $\kappa$ is the gradient coefficient. Note that, for simplicity, we assume fast diffusion and only keep the reaction term.

Under galvanostatic condition, we constraint the total reaction rate by an imposed rate $I_{\text{app}}$:

$$
\begin{equation}
\int_{\Omega} R(\eta) \, d\Omega = I_{\text{app}}.
\end{equation}
$$


### Initial Condition

The initial condition is similar to `allen_cahn` application, but rescaled between c1 and c2.

### Boundary Condition

Natural boundary condition is used here.


## Weak form
Here, we chose $\mu$ as a separate variable. We multiply the PDE by a sufficiently smooth test function $w(\vec{x})$ and integrate over the entire domain $\Omega$:

$$
\begin{equation}
\int_{\Omega} w ~ \frac{\partial c}{\partial t} \, d\Omega = \int_{\Omega} w ~ R(\mu+V) \, d\Omega
\end{equation}
$$

$$
\begin{equation}
\int_{\Omega} w ~ \mu \, d\Omega = \int_{\Omega} w ~ \left( f_{,c}(c) - \kappa \nabla^2 c \right) \, d\Omega
\end{equation}
$$

Using Integration by Parts and the Divergence Theorem:

$$
\begin{equation}
\int_{\Omega} w ~ \mu \, d\Omega = \int_{\Omega} w ~ f_{,c}(c) + \nabla w \cdot (\kappa\nabla c) \, d\Omega
\end{equation}
$$

## Time discretization: Forward Euler (explicit)
We discretize the time derivative by

$$
\begin{equation}
\frac{\partial c}{\partial t} = \frac{c^{n} - c^{n-1}}{\Delta t}
\end{equation}
$$

and use $c^{n-1}$ for the right hand side of the equation.

The weak form becomes

$$
\begin{equation}
\int_{\Omega} w ~ c^{n} \, d\Omega = \int_{\Omega} w \left(c^{n-1} + \Delta t R(\eta^{n-1})\right) \, d\Omega
\end{equation}
$$

$$
\begin{equation}
\int_{\Omega} w ~ \mu^{n} \, d\Omega = \int_{\Omega} w \left(f_{,c}(c^{n})\right) + \nabla w \cdot (\kappa \nabla c^{n}) \, d\Omega
\end{equation}
$$

We have

$$
\begin{equation}
r_{c}= c^{n-1} + \Delta t ~ R(\mu^{n-1}+V^{n-1})
\end{equation}
$$

$$
\begin{equation}
r_{\mu}= f_{,c}(c^{n})
\end{equation}
$$

$$
\begin{equation}
r_{\eta x} = \kappa \nabla c^{n}
\end{equation}
$$

## Galvanostatic condition
The applied voltage $V^n$ need to be solved from the galvanostatic condition (Eq. 2). Define the residual

$$
r(V^n) = \int_{\Omega} R(\eta^n(\mu^n, V^n)) \, d\Omega - I_{\text{app}}.
$$

This is a nonlinear equation of $V^n$, which can be solved using Newton's method (start from an initial guess $V_0$)

$$
V^n = V_0 - \frac{r(V_0)}{r'(V_0)},
$$

We then use the updated $V^n$ as the new initial guess and repeat until convergence.

## Special function form
We choose the regular solution free energy density

$$
f(c) = m c(1-c) + c \ln c + (1-c) \ln (1-c)
$$

Its derivative is

$$
f_{,c}(c) = m (1-2c) + \ln \frac{c}{1-c}
$$

The reaction kinetics

$$
R(\eta) = 2 \sinh \frac{\eta}{2}
$$

Its derivative is

$$
R_{,\eta}(\eta) = \cosh \frac{\eta}{2}
$$
