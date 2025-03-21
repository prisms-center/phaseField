## PRISMS-PF: Allen-Cahn Dynamics
Consider a free energy expression of the form:

$$
\begin{align}
\Pi(\eta, \nabla \eta) = \int_{\Omega} f(\eta) + \frac{\kappa}{2} \nabla \eta \cdot \nabla \eta ~dV 
\end{align}
$$

where $\eta$ is the structural order parameter, and $\kappa$ is the gradient length scale parameter.

### Variational treatment
Considering variations on the primal field $\eta$ of the from $\eta+\epsilon w$, we have

$$
\begin{align}
\delta \Pi &=  \left. \frac{d}{d\epsilon} \int_{\Omega}  f(\eta+\epsilon w) +  \frac{\kappa}{2} \nabla  (\eta+\epsilon w)  \cdot  ~\nabla  (\eta+\epsilon w)   ~dV \right\vert_{\epsilon=0} \\
&=  \int_{\Omega}   w f_{,\eta} +   \kappa \nabla w \nabla  \eta    ~dV \\
&=  \int_{\Omega}   w \left( f_{,\eta} -  \kappa \Delta \eta \right)  ~dV  +   \int_{\partial \Omega}   w \kappa \nabla \eta \cdot n   ~dS
\end{align}
$$

Assuming $\kappa \nabla \eta \cdot n = 0$, and using standard variational arguments on the equation $\delta \Pi =0$, we have the expression for chemical potential as

$$
\begin{align}
\mu  = f_{,\eta} -  \kappa \Delta \eta
\end{align}
$$

## Kinetics
Now the parabolic PDE for Allen-Cahn dynamics is given by:

$$
\begin{align}
\frac{\partial \eta}{\partial t} = -M(f_{,\eta} - \kappa \Delta \eta)
\end{align}
$$

where $M$ is the constant mobility.

## Time discretization

Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:

$$
\begin{align}
\eta^{n+1} = \eta^{n} - \Delta t M(f_{,\eta}^{n} - \kappa \Delta \eta^{n})
\end{align}
$$

## Weak formulation

In the weak formulation, considering an arbitrary variation $w$, the above equation can be expressed as a residual equation:

$$
\begin{align}
\int_{\Omega} w \eta^{n+1} ~dV &= \int_{\Omega} w \eta^{n} - w \Delta t M(f_{,\eta}^{n} - \kappa \Delta \eta^{n}) ~dV \\
&= \int_{\Omega} w \left( \eta^{n} - \Delta t M f_{,\eta}^{n} \right) + \nabla w (-\Delta t M \kappa) \cdot (\nabla \eta^{n}) ~dV \quad [\kappa \nabla \eta \cdot n = 0 \quad \text{on} \quad \partial \Omega]
\end{align}
$$

$$
\begin{align}
r_{\eta}= \eta^{n} - \Delta t M f_{,\eta}^{n} 
\end{align}
$$

$$
\begin{align}
r_{\eta x} = -\Delta t M \kappa\nabla \eta^{n}
\end{align}
$$

The above values of $r_{\eta}$ and $r_{\eta x}$ are used to define the residuals in `applications/allenCahn/equations.cc`.