# PRISMS PhaseField: Cahn-Hilliard Dynamics (Mixed-Formulation)
Consider a free energy expression of the form:

\begin{equation}
  \Pi(c, \grad  c) = \int_{\Omega}    f( c ) + \frac{\kappa}{2} \grad  c  \cdot \grad  c    ~dV 
\end{equation}
where $c$ is the composition, and $\kappa$ is the gradient length scale parameter.
	
\section{Variational treatment}
Considering variations on the primal field $c$ of the from $c+\epsilon w$, we have
\begin{align}
\delta \Pi &=  \left. \frac{d}{d\epsilon} \int_{\Omega}  f(c+\epsilon w) +  \frac{\kappa}{2} \grad  (c+\epsilon w)  \cdot  ~\grad  (c+\epsilon w)   ~dV \right\vert_{\epsilon=0} \\
&=  \int_{\Omega}   w f_{,c} +   \kappa \grad w \grad  c    ~dV \\
&=  \int_{\Omega}   w \left( f_{,c} -  \kappa \Delta c \right)  ~dV  +   \int_{\partial \Omega}   w \kappa \grad c \cdot n   ~dS
\end{align}
Assuming $\kappa \grad c \cdot n = 0$, and using standard variational arguments on the equation $\delta \Pi =0$ we have the expression for chemical potential as
\begin{equation}
  \mu  = f_{,c} -  \kappa \Delta c
\end{equation}

\section{Kinetics}
Now the Parabolic PDE for Cahn-Hilliard dynamics is given by:
\begin{align}
  \frac{\partial c}{\partial t} &= -~\grad \cdot (-M\grad \mu)\\
  &=-M~\grad \cdot (-\grad (f_{,c} -  \kappa \Delta c)) 
\end{align}
where $M$ is the constant mobility. This equation can be split into two equations as follow:
\begin{align}
  \mu &= f_{,c} -  \kappa \Delta c \\
  \frac{\partial c}{\partial t} &= M~\grad \cdot (\grad \mu)
\end{align}
\section{Time discretization}
Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:
\begin{align}
  \mu^{n+1} &= f_{,c}^{n} -  \kappa \Delta c^{n} \\
 c^{n+1} &= c^{n} + \Delta t M~\grad \cdot (\grad \mu^{n})
\end{align}

\section{Weak formulation}
In the weak formulation, considering an arbitrary variation $w$, the above equations can be expressed as residual equations representing a mixed (split) formulation:
\begin{align}
  \int_{\Omega}   w  \mu^{n+1}  ~dV &= \int_{\Omega}  w  f_{,c}^{n} - w \kappa \Delta c^{n}~dV\\
  &=\int_{\Omega}  w  \underbrace{f_{,c}^{n}}_{r_{mu}} + \grad w \cdot \underbrace{\kappa \grad c^{n}}_{r_{mux}} ~dV 
\end{align}
and 
\begin{align}
\int_{\Omega}   w c^{n+1} ~dV&= \int_{\Omega}   w c^{n} + w \Delta t M~\grad \cdot (\grad \mu^{n}) ~dV \\
&= \int_{\Omega}   w \underbrace{ c^{n}}_{r_{c}} + \grad w \underbrace{ (-\Delta t M)~ \cdot (\grad \mu^{n})}_{r_{c x}} ~dV \quad \text{[neglecting boundary flux]} 
\end{align}
\vskip 0.25in
The above values of $r_{mu}$, $r_{mux}$, $r_{c}$ and $r_{cx}$ are used to define the residuals in the following parameters file: \\
\textit{applications/cahnHilliard/parameters.h}
