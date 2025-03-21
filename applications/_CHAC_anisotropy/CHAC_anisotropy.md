# PRISMS-PF: CHAC Anisotropy (with Coupled CH-AC Dynamics)
Consider a free energy expression of the form:

$$
\begin{equation}
  \Pi(c, \eta, \nabla  \eta) = \int_{\Omega}    \left( f_{\alpha}(1-H) + f_{\beta}H \right)  + \frac{1}{2} | \gamma( \mathbf{n} ) \nabla  \eta |^2   ~dV 
\end{equation}
$$

where $f_{\alpha}$ and $f_{\beta}$ are the free energy densities corresponding to $\alpha$ and $\beta$ phases, respectively, and are functions of composition $c$. $H$ is a function of the structural order parameter $\eta$.  The interface normal vector $\mathbf{n}$ is given by 

$$
\begin{equation}
\mathbf{n} = \frac{\nabla \eta}{|\nabla \eta|}
\end{equation}
$$

for $\nabla \eta \ne \mathbf{0}$, and $\mathbf{n} = \mathbf{0}$ when $\nabla \eta = \mathbf{0}$.

## Variational Treatment

Following standard variational arguments (see Cahn-Hilliard formulation), we obtain the chemical potentials:

$$
\begin{align}
  \mu_{c}  &= (f_{\alpha,c}(1-H)+f_{\beta,c}H)  \\
  \mu_{\eta}  &= (f_{\beta,c}-f_{\alpha,c})H_{,\eta} - \nabla \cdot \mathbf{m}
\end{align}
$$

The components of the anisotropic gradient $\mathbf{m}$ are given by

$$
\begin{equation}
m_i = \gamma(\mathbf{n}) \left( \nabla \eta + |\nabla \eta| (\delta_{ij}-n_i n_j) \frac{\partial \gamma (\mathbf{n})}{n_j} \right),
\end{equation}
$$

where $\delta_{ij}$ is the Kronecker delta.

## Kinetics

Now the PDE for Cahn-Hilliard dynamics is given by:

$$
\begin{align}
  \frac{\partial c}{\partial t} &= -~\nabla \cdot (-M_c\nabla \mu_c)\\
  &=M_c~\nabla \cdot (\nabla (f_{\alpha,c}(1-H)+f_{\beta,c}H)) 
  \end{align}
$$

and the PDE for Allen-Cahn dynamics is given by:

$$
\begin{align}
  \frac{\partial \eta}{\partial t} &= -M_\eta \mu_\eta \\
  &=-M_\eta ~ ((f_{\beta,c}-f_{\alpha,c})H_{,\eta} - \nabla \cdot \mathbf{m}) 
\end{align}
$$

where $M_c$ and $M_\eta$ are the constant mobilities. 

## Time discretization

Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:

$$
\begin{align}
 \eta^{n+1} &= \eta^{n}  - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n -  \nabla \cdot \mathbf{m}^n) \\
 c^{n+1} &= c^{n}  + \Delta t M_{\eta}~\nabla \cdot (\nabla (f_{\alpha,c}^n(1-H^{n})+f_{\beta,c}^n H^{n}))
\end{align}
$$

## Weak formulation
In the weak formulation, considering an arbitrary variation $w$, the above equations can be expressed as residual equations:

$$
\begin{align}
\int_{\Omega}   w  \eta^{n+1}  ~dV &= \int_{\Omega}   w \eta^{n} -   w    \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n - \kappa \Delta \eta^n)  ~dV
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}  w  \left( \eta^{n} - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n) \right)+ \nabla w \cdot (- \Delta t M_{\eta})   \mathbf{m}^n ~dV 
\end{align}
$$

$$
\begin{align}
r_{\eta} &= \eta^{n} - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n)
\end{align}
$$

$$
\begin{align}
r_{\eta x} &= (- \Delta t M_{\eta})   \mathbf{m}^n
\end{align}
$$

and 

$$
\begin{align}
\int_{\Omega}   w  c^{n+1}  ~dV &= \int_{\Omega}   w c^{n} + w    \Delta t M_{c}~ \nabla \cdot (\nabla (f_{\alpha,c}^n(1-H^{n})+f_{\beta,c}^n H^{n}))  ~dV\\
&= \int_{\Omega}   w c^{n} +  \nabla w   (-\Delta t M_{c})~ [~(f_{\alpha,cc}^n(1-H^{n})+f_{\beta,cc}^n H^{n}) \nabla c + ~((f_{\beta,c}^n-f_{\alpha,c}^n)H^{n}_{,\eta} \nabla \eta) ] ~dV
\end{align}
$$

$$
\begin{align}
r_c &= c^{n}
\end{align}
$$

$$
\begin{align}
r_{cx} &= (-\Delta t M_{c})~ [~(f_{\alpha,cc}^n(1-H^{n})+f_{\beta,cc}^n H^{n}) \nabla c + ~((f_{\beta,c}^n-f_{\alpha,c}^n)H^{n}_{,\eta} \nabla \eta) ] 
\end{align}
$$

The above values of $r_{\eta}$, $r_{\eta x}$, $r_{c}$ and $r_{cx}$ are used to define the residuals in the following equations file: 
`applications/CHAC\_anisotropy/equations.h`
