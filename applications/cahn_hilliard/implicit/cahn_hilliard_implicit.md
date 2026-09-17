# PRISMS PhaseField: Cahn-Hilliard Dynamics (Mixed-Formulation)
Consider a free energy expression of the form:

$$
\begin{equation}
  \Pi(c, \nabla  c) = \int_{\Omega}    f( c ) + \frac{\kappa}{2} \nabla  c  \cdot \nabla  c    ~dV
\end{equation}
$$

where $c$ is the composition, and $\kappa$ is the gradient length scale parameter.

## Variational treatment
Considering variations on the primal field $c$ of the from $c+\epsilon w$, we have

$$
\begin{align}
\delta \Pi &=  \left. \frac{d}{d\epsilon} \int_{\Omega}  f(c+\epsilon w) +  \frac{\kappa}{2} \nabla  (c+\epsilon w)  \cdot  ~\nabla  (c+\epsilon w)   ~dV \right\vert_{\epsilon=0}
\end{align}
$$

$$
\begin{align}
&=  \int_{\Omega}   w f_{,c} +   \kappa \nabla w \nabla  c    ~dV
\end{align}
$$

$$
\begin{align}
&=  \int_{\Omega}   w \left( f_{,c} -  \kappa \nabla^2 c \right)  ~dV  +   \int_{\partial \Omega}   w \kappa \nabla c \cdot n   ~dS
\end{align}
$$

Assuming $\kappa \nabla c \cdot n = 0$, and using standard variational arguments on the equation $\delta \Pi =0$ we have the expression for chemical potential as

$$
\begin{equation}
  \mu  = f_{,c} -  \kappa \nabla^2 c
\end{equation}
$$

## Kinetics
Now the Parabolic PDE for Cahn-Hilliard dynamics is given by:

$$
\begin{align}
  \frac{\partial c}{\partial t} &= -~\nabla \cdot (-M\nabla \mu)
\end{align}
$$

$$
\begin{align}
  &=-M~\nabla \cdot (-\nabla (f_{,c} -  \kappa \nabla^2 c))
\end{align}
$$

where $M$ is the constant mobility. This equation can be split into two equations as follow:

$$
\begin{align}
  \mu &= f_{,c} -  \kappa \nabla^2 c
\end{align}
$$

$$
\begin{align}
  \frac{\partial c}{\partial t} &= \nabla \cdot (M\nabla \mu)
\end{align}
$$

## Time discretization

Considering backward Euler explicit time stepping, we have the time discretized kinetics equation:

$$
\begin{align}
  \frac{c^{n} - c^{n-1}}{\Delta t} &= \nabla \cdot (M\nabla \mu^{n-1})
\end{align}
$$

and

$$
\begin{align}
  \mu^{n} &= f_{,c}^{n} -  \kappa \nabla^2 c^{n}
\end{align}
$$

## Newton Form

To solve this system of equations using Newton's method we construct residual expressions by rearranging equations 10 and 11

$$
\begin{align}
  0 = R_c &= c^{n-1} - c^{n} + \Delta t~\nabla \cdot (M\nabla \mu^{n})  \\
  &= c^{n-1} - c^{n} + \nabla \cdot (\Delta tM\nabla \mu^{n})
\end{align}
$$

and

$$
\begin{align}
  0 = R_{\mu} &= \mu^{n} - f_{,c}^{n} + \kappa \nabla^2 c^{n} \\
  &= \mu^{n} - f_{,c}^{n} + \nabla \cdot \kappa \nabla c^{n}
\end{align}
$$

We also need the linearization of the residual (the Jacobian) evaluated at a test change.

$$
\begin{align}
  \frac{\partial R_c}{\partial c^{n}}\bigg\vert_{\Delta c} 
  &= -\Delta c
\\
  \frac{\partial R_c}{\partial \mu^{n}}\bigg\vert_{\Delta \mu} 
  &= \nabla \cdot (\Delta tM\nabla (\Delta \mu))
\end{align}
$$

$$
\begin{align}
  \frac{\partial R_\mu}{\partial c^{n}}\bigg\vert_{\Delta c} 
  &= - f_{,cc}^{n}\Delta c + \nabla \cdot \kappa \nabla (\Delta c)
\\
  \frac{\partial R_\mu}{\partial \mu^{n}}\bigg\vert_{\Delta \mu} 
  &= \Delta \mu
\end{align}
$$

## Weak Formulation

For c,

$$
\begin{align}
\int_{\Omega} w (-( - \Delta c + \nabla \cdot (\Delta tM\nabla (\Delta \mu)))) ~dV
&=
\int_{\Omega} w (c^{n-1} - c^{n} + \nabla \cdot (\Delta tM\nabla \mu^{n})) ~dV
\\

\int_{\Omega} w ( \Delta c) + \nabla w\cdot (\Delta tM\nabla (\Delta \mu)) ~dV
&=
\int_{\Omega} w (c^{n-1} - c^{n}) + \nabla w \cdot (-\Delta tM\nabla \mu^{n}) ~dV
\end{align}
$$

with 

$$
\begin{align}
LHS_{c}   &= \Delta c \\
LHS_{cx}  &= \Delta tM\nabla (\Delta \mu) \\
RHS_{c}   &= c^{n-1} - c^{n} \\
RHS_{c x} &= -\Delta tM\nabla \mu^{n} \\
\end{align}
$$

and for $\mu$,

$$
\begin{align}
\int_{\Omega} w (-(\Delta \mu- f_{,cc}^{n}\Delta c + \nabla \cdot \kappa \nabla (\Delta c))) ~dV
&=
\int_{\Omega} w (\mu^{n} - f_{,c}^{n} + \nabla \cdot \kappa \nabla c^{n}) ~dV
\\
\int_{\Omega} w (-\Delta \mu + f_{,cc}^{n}\Delta c) + \nabla w \cdot (\kappa \nabla (\Delta c)) ~dV
&=
\int_{\Omega} w (\mu^{n} - f_{,c}^{n}) +\nabla w \cdot (-\kappa \nabla c^{n}) ~dV
\end{align}
$$

with 

$$
\begin{align}
LHS_{\mu }   &= -\Delta \mu + f_{,cc}^{n}\Delta c \\
LHS_{\mu x}  &= \kappa \nabla (\Delta c) \\
RHS_{\mu }   &= \mu^{n} - f_{,c}^{n} \\
RHS_{\mu  x} &= -\kappa \nabla c^{n} \\
\end{align}
$$


The above expressions define the code written in:
`applications/cahn_hilliard/implicit/custom_pde.h`
