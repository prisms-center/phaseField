# PRISMS-PF: Eshelby Inclusion

This example application implements a simple 3D calculation of the displacement field near a homogeneous inclusion.

Consider a strain energy expression of the form:

$$
\begin{equation}
  \Pi(\varepsilon) = \int_{\Omega}  \frac{1}{2} (\epsilon-\epsilon^0):C: (\epsilon-\epsilon^0)   ~dV  - \int_{\partial \Omega}   u \cdot t  ~dS
\end{equation}
$$

where $\varepsilon$ is the infinitesimal strain tensor, $C_{ijkl}=\lambda \delta_{ij} \delta_{kl}+\mu ( \delta_{ik} \delta_{jl}+ \delta_{il} \delta_{jk} )$  is the fourth order elasticity tensor, ($\lambda$, $\mu$) are the Lame parameters, $u$ is the displacement, and $t=\sigma \cdot n$ is the surface traction. Here, body forces are assumed to be zero.

## Governing equations
Considering variations on the displacement $u$ of the from $u+\alpha w$, we have

$$
\begin{align}
\delta \Pi &=  \left. \frac{d}{d\alpha} \left( \int_{\Omega}   \frac{1}{2} [\epsilon(u+\alpha w) - \epsilon^0] : C : [\epsilon(u+\alpha w) - \epsilon^0] ~dV - \int_{\partial \Omega} (u+\alpha w) \cdot t ~dS \right)\right\vert_{\alpha=0} \\
&=  \int_{\Omega}   \nabla w : C :  (\epsilon-\epsilon^0)  ~dV -  \int_{\partial \Omega}   w \cdot t~dS\\
&=  \int_{\Omega}   \nabla w : \sigma  ~dV -  \int_{\partial \Omega}   w\cdot t~dS
\end{align}
$$

where $\sigma = C : (\epsilon-\epsilon^0)$ is the stress tensor, $\epsilon^0$ is the misfit strain (eigenstrain), and $t=\sigma \cdot n$ is the surface traction. In this case, we assume that the diagonal elements of $\epsilon^0$ take the form:

$$
\begin{equation}
\epsilon^0_{ii} = m \left(\frac{1}{2} + \frac{1}{2} \tanh\left(\frac{r-a}{\ell}\right) \right)
\end{equation}
$$

where $m$ is the magnitide of the misfit strain inside the inclusion, $\ell$ determines the thickness of the "interface" between the inclusion and the matrix, $r$ is the distance from the origin, and $a$ is the radius of the inclusion. The off-diagonal elements of $\epsilon^0$ are zero.

The minimization of the variation, $\delta \Pi=0$, gives the weak formulation of the governing equation of mechanics:

$$
\begin{align}
\int_{\Omega}   \nabla w : \sigma  ~dV -  \int_{\partial \Omega}   w \cdot t  ~dS = 0
\end{align}
$$

If surface tractions are zero:

$$
\begin{align}
R &=  \int_{\Omega}   \nabla w :  \sigma ~dV = 0
\end{align}
$$

## Residual expressions
In PRISMS-PF, two sets of residuals are required for elliptic PDEs (such as this one), one for the left-hand side of the equation (LHS) and one for the right-hand side of the equation (RHS). We solve $R=0$ by casting this in a form that can be solved as a matrix inversion problem. This will involve a brief detour into the discretized form of the equation. First we derive an expression for the solution, given an initial guess, $u_0$:

$$
\begin{gather}
0 = R(u) = R(u_0 + \Delta u)
\end{gather}
$$

where $\Delta u = u - u_0$. Then, we expand $R(u)$ around $u_0$:
$$
\begin{equation}
R(u) = R(u_0 + \Delta u) \approx R(u_0) + \int_\Omega \left.\frac{\delta R}{\delta u}\right|_{u_0} \cdot \Delta u ~dV
\end{equation}
$$

The variation of $R$ is
$$
\begin{align}
\delta R &= \int_{\Omega} \nabla w : C : \nabla \delta u ~dV\\
         &= - \int_{\Omega} \nabla \cdot (\nabla w : C) \cdot \delta u ~dV
\end{align}
$$

Therefore, the variational derivative of $R$ is
$$
\begin{equation}
\frac{\delta R}{\delta u} = -\nabla \cdot (\nabla w : C)
\end{equation}
$$

The residual equation can be written as:
$$
\begin{align}
R(u_0) &= \int_\Omega (\nabla \cdot (\nabla w : C)) \cdot \Delta u ~dV\\
&= -\int_{\Omega} \nabla w : C : \nabla (\Delta u) dV
\end{align}
$$

Thus, the full equation relating $u_0$ and $\Delta u$ is:

$$
\begin{align}
\int_{\Omega} \nabla w : C : \nabla (\Delta u) dV = -\int_{\Omega}   \nabla w : \sigma ~dV
\end{align}
$$

$$
\begin{align}
r_{ux}^{LHS} &= C : \nabla (\Delta u)
\end{align}
$$

$$
\begin{align}
r_{ux} &= -\sigma
\end{align}
$$

The above values of $r_{ux}^{LHS}$ and $r_{ux}$ are used to define the residuals in the following input file:
`applications/eschelby_inclusion/equations.cc`
