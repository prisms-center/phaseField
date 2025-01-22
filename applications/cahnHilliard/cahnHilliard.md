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
&=  \int_{\Omega}   w \left( f_{,c} -  \kappa \Delta c \right)  ~dV  +   \int_{\partial \Omega}   w \kappa \nabla c \cdot n   ~dS
\end{align}
$$

Assuming $\kappa \nabla c \cdot n = 0$, and using standard variational arguments on the equation $\delta \Pi =0$ we have the expression for chemical potential as

$$
\begin{equation}
  \mu  = f_{,c} -  \kappa \Delta c
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
  &=-M~\nabla \cdot (-\nabla (f_{,c} -  \kappa \Delta c)) 
\end{align}
$$

where $M$ is the constant mobility. This equation can be split into two equations as follow:

$$
\begin{align}
  \mu &= f_{,c} -  \kappa \Delta c
\end{align}
$$

$$
\begin{align}
  \frac{\partial c}{\partial t} &= M~\nabla \cdot (\nabla \mu)
\end{align}
$$

## Time discretization

Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:

$$
\begin{align}
  \mu^{n+1} &= f_{,c}^{n} -  \kappa \Delta c^{n} 
\end{align}
$$

$$
\begin{align}
 c^{n+1} &= c^{n} + \Delta t M~\nabla \cdot (\nabla \mu^{n})
\end{align}
$$

## Weak formulation
In the weak formulation, considering an arbitrary variation $w$, the above equations can be expressed as residual equations representing a mixed (split) formulation:

$$
\begin{align}
  \int_{\Omega}   w  \mu^{n+1}  ~dV &= \int_{\Omega}  w  f_{,c}^{n} - w \kappa \Delta c^{n}~dV
\end{align}
$$

$$
\begin{align}
 &=\int_{\Omega}  w  f_{,c}^{n} + \nabla w \cdot \kappa \nabla c^{n} ~dV  
\end{align}
$$

$$
\begin{align}
r_{mu} &= f_{,c}^{n}
\end{align}
$$

$$
\begin{align}
r_{mux} &= \kappa \nabla c^{n}
\end{align}
$$

and 

$$
\begin{align}
\int_{\Omega}   w c^{n+1} ~dV&= \int_{\Omega}   w c^{n} + w \Delta t M~\nabla \cdot (\nabla \mu^{n}) ~dV 
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}   w c^{n} + \nabla w  (-\Delta t M)~ \cdot (\nabla \mu^{n}) ~dV \quad \text{[neglecting boundary flux]} 
\end{align}
$$

$$
\begin{align}
r_{c} &= c^{n}
\end{align}
$$

$$
\begin{align}
r_{c x} &= (-\Delta t M)~ \cdot (\nabla \mu^{n})
\end{align}
$$

The above values of $r_{mu}$, $r_{mux}$, $r_{c}$ and $r_{cx}$ are used to define the residuals in the following parameters file: 
`applications/cahnHilliard/parameters.h`
