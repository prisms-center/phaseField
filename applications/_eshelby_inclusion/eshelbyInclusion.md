# PRISMS-PF Application Formulation: eshelbyInclusion

This example application implements a simple 3D calculation of the displacement field near a homogenous inclusion. 

Consider a strain energy expression of the form:

$$
\begin{equation}
  \Pi(\varepsilon) = \int_{\Omega}  \frac{1}{2} (\epsilon-\epsilon^0):C: (\epsilon-\epsilon^0)   ~dV  - \int_{\partial \Omega}   u \cdot t  ~dS
\end{equation}
$$

where $\varepsilon$ is the infinitesimal strain tensor, $C_{ijkl}=\lambda \delta_{ij} \delta_{kl}+\mu ( \delta_{ik} \delta_{jl}+ \delta_{il} \delta_{jk} )$  is the fourth order elasticity tensor, ($\lambda$, $\mu$) are the Lame parameters, $u$ is the displacement, and $t=\sigma \cdot n$ is the surface traction.

## Governing equations
Considering variations on the displacement $u$ of the from $u+\alpha w$, we have

$$
\begin{align}
\delta \Pi &=  \left. \frac{d}{d\alpha} \left( \int_{\Omega}   \frac{1}{2} [\epsilon(u+\alpha w) - \epsilon^0] : C : [\epsilon(u+\alpha w) - \epsilon^0] ~dV - \int_{\partial \Omega} u \cdot t ~dS \right)\right\vert_{\alpha=0} \\
&=  \int_{\Omega}   \nabla w : C :  (\epsilon-\epsilon^0)  ~dV -  \int_{\partial \Omega}   w \cdot t~dS\\
&=  \int_{\Omega}   \nabla w : \sigma  ~dV -  \int_{\partial \Omega}   w\cdot t~dS
\end{align}
$$

where $\sigma = C : (\epsilon-\epsilon^0)$ is the stress tensor, $\epsilon^0$ is the misfit strain (eigenstrain), and $t=\sigma \cdot n$ is the surface traction. In this case, we assume that the diagonal elements of $\epsilon^0$ take the form:

$$
\begin{equation}
\epsilon^0_{ii} = m \left(\frac{1}{2} + \frac{1}{2} \tanh(l(r-a)) \right)
\end{equation}
$$

where $m$ is the magnitide of the misfit strain inside the inclusion, $l$ determines the thickness of the ``interface'' between the inclusion and the matrix, $r$ is the distance from the origin, and $a$ is the radius of the inclusion. The off-diagonal elements of $\epsilon^0$ are zero.

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

where $\Delta u = u - u_0$. Then, applying the discretization that $u = \sum_i w^i U^i$, we can write the following linearization:

$$
\begin{equation}
\frac{\delta R(u)}{\delta u} \Delta U = -R(u_0) 
\end{equation}
$$

The discretized form of this equation can be written as a matrix inversion problem. However, in PRISMS-PF, we only care about the product $\frac{\delta R(u)}{\delta u} \Delta U$. Taking the variational derivative of $R(u)$ yields:

$$
\begin{align}
\frac{\delta R(u)}{\delta u} &= \frac{d}{d\alpha} \int_{\Omega}   \nabla w :C: \left[ \epsilon (u+\alpha w) - \epsilon^0 \right] ~dV  \bigg{|}_{\alpha=0} 
\end{align}
$$

$$
\begin{align}
&=  \int_{\Omega}   \nabla w :C: \frac{1}{2}\frac{d}{d\alpha}\left[ \nabla(u+\alpha w) + \nabla(u+\alpha w)^T  - \epsilon^0\right] ~dV \bigg{|}_{\alpha=0}
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}   \nabla w :C: \frac{d}{d\alpha} \left[ \nabla(u+\alpha w) - \epsilon^0 \right]  ~dV \bigg{|}_{\alpha=0} \quad (due ~to ~the ~symmetry ~of ~C) 
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}   \nabla w :C: \nabla w  ~dV 
\end{align}
$$

In its discretized form $\frac{\delta R(u)}{\delta u} \Delta U$ is:

$$
\begin{equation}
\frac{\delta R(u)}{\delta u} \Delta U = \sum_i \sum_j \int_{\Omega} \nabla N^i : C : \nabla N^j dV ~\Delta U^j
\end{equation}
$$

Moving back to the non-discretized form yields:

$$
\begin{equation}
\frac{\delta R(u)}{\delta u} \Delta U = \int_{\Omega} \nabla w : C : \nabla (\Delta u) dV
\end{equation}
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
r_{ux} &= \sigma
\end{align}
$$

The above values of $r_{ux}^{LHS}$ and $r_{ux}$ are used to define the residuals in the following input file:
`applications/eschelbyInclusion/equations.h`
