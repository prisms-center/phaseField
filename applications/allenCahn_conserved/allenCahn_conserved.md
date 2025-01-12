# PRISMS PhaseField: Globally Conserved Allen-Cahn Dynamics}}

This application performs a  phase field simulation of Allen-Cahn dynamics subject to global (as opposed to local) conservation. Global conservation implies that $\int_\Omega \eta \,dV$, where $\Omega$ is the volume of the system, is constant in time.

Note: This application needs to run with **uniform mesh** in order for the calculation of the chemical potential to be accurate.

Consider a free energy expression of the form:

$$
\begin{equation}
  \Pi(\eta, \nabla  \eta) = \int_{\Omega}    f( \eta ) + \frac{\kappa}{2} \nabla  \eta  \cdot \nabla  \eta    ~dV, 
\end{equation}
$$

where $\eta$ is the structural order parameter, and $\kappa$ is the gradient length scale parameter.
	
## Variational treatment

Considering variations on the primal field $\eta$ of the from $\eta+\epsilon w$, we have

$$
\begin{align}
\delta \Pi &=  \left. \frac{d}{d\epsilon} \int_{\Omega}  f(\eta+\epsilon w) +  \frac{\kappa}{2} \nabla  (\eta+\epsilon w)  \cdot  ~\nabla  (\eta+\epsilon w)   ~dV \right\vert_{\epsilon=0}
\end{align}
$$

$$
\begin{align}
&=  \int_{\Omega}   w f_{,\eta} +   \kappa \nabla w \nabla  \eta    ~dV
\end{align}
$$

$$
\begin{align}
&=  \int_{\Omega}   w \left( f_{,\eta} -  \kappa \Delta \eta \right)  ~dV  +   \int_{\partial \Omega}   w \kappa \nabla \eta \cdot n   ~dS,
\end{align}
$$

where $f_{,\eta} = \partial f/\partial \eta$.

Assuming $\kappa \nabla \eta \cdot n = 0$, and using standard variational arguments on the equation $\delta \Pi =0$ we have the expression for chemical potential as

$$
\begin{equation}
  \mu  = f_{,\eta} -  \kappa \Delta \eta.
\end{equation}
$$

## Kinetics
The Parabolic PDE for Allen-Cahn dynamics is given by:

$$
\begin{align}
\frac{\partial \eta}{\partial t} &= -M\mu,
\end{align}
$$

where $M$ is the constant mobility. However, the previous equation does not ensure global conservation of $\eta$. In order to achieve global conservation we can add a spatially-uniform term, $A$, to the RHS that effectively offsets the total change in $\eta$ such that $\int_\Omega \partial \eta / \partial t\,dV =0$:

$$
\begin{align}
\frac{\partial \eta}{\partial t} &= -M\mu+A.
\end{align}
$$

Applying the global conservation constraint to previous equation, we obtain

$$
\begin{align}
\int_\Omega \frac{\partial \eta}{\partial t}\,dV &= -\int_\Omega (M\mu-A) \,dV = 0.
\end{align}
$$

Since $A$ is spatially uniform and $M$ is a constant, $A$ given by

$$
\begin{align}
A &= \frac{M}{V}\int_\Omega \mu \,dV.
\end{align}
$$

Subsituting $A$ from

$$
\begin{align}
A &= \frac{M}{V}\int_\Omega \mu \,dV.
\end{align}
$$

into 

$$
\begin{align}
\frac{\partial \eta}{\partial t} &= -M\mu+A.
\end{align}
$$

we get

$$
\begin{align}
\frac{\partial \eta}{\partial t} &= -M(\mu-\bar{\mu}),
\end{align}
$$

where $\bar{\mu} = (1/V)\int_\Omega \mu \,dV$.

## Time discretization
Considering forward Euler explicit time stepping, we have the time discretized kinetics equations:

$$
\begin{align}
\eta^{n+1} &= \eta^{n} - \Delta t M(\mu^n-\bar{\mu}^n)
\end{align}
$$

and

$$
\begin{align}
\mu^{n+1} &= f_{,\eta}^n -  \kappa \Delta \eta^n.
\end{align}
$$

## Weak formulation
In the weak formulation, considering an arbitrary variation $w$, the above equation can be expressed as a residual equations:

$$
\begin{align}
\int_{\Omega}   w \eta^{n+1} ~dV&= \int_{\Omega}  w [ \eta^{n} - \Delta t M(\mu^n-\bar{\mu}^n) ]~dV 
\end{align}
$$

$$
\begin{align}
r_{\eta} &= \eta^{n} - \Delta t M(\mu^n-\bar{\mu}^n)
\end{align}
$$

and

$$
\begin{align}
\int_{\Omega}   w \mu^{n+1} ~dV&= \int_{\Omega} [w f_{,\eta}^n - w \kappa \Delta \eta^{n}]~dV 
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}   w (  f_{,\eta}^{n} )  + \nabla w \cdot (\kappa \nabla \eta^{n}) ~dV, 
\end{align}
$$

$$
\begin{align}
r_{\mu x} &=  (\kappa \nabla \eta^{n})
\end{align}
$$

$$
\begin{align}
r_{\mu} &=  f_{,\eta}^{n} 
\end{align}
$$

where the reference chemical potential $\bar{\mu}^n$ for time step $n$ in

$$
\begin{align}
\int_{\Omega}   w \eta^{n+1} ~dV&= \int_{\Omega}  w [ \eta^{n} - \Delta t M(\mu^n-\bar{\mu}^n) ]~dV 
\end{align}
$$

is calculated as

$$
\begin{align}
\bar{\mu}^n& = \frac{1}{V}\int_\Omega \mu^n ~dV.
\end{align}
$$

The above values of  $r_{\eta}$, $r_{\mu}$, and $r_{\mu x}$ are used to define the residuals in the following parameters file:
`\textit{applications/allenCahn\_conserved/equations.cc`
