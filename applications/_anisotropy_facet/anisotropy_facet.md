# PRISMS PhaseField: Faceted Anisotropy (with Coupled CH-AC Dynamics)

This application is essentially a specialization of the CHAC\_anisotropyRegarized application with a particular choice of interfacial energy anisotropy $\gamma(\mathbf{n})$.  In this document, we repeat the formulation for that model for completeness and then describe the anisotropy used in this application.  Consider a free energy expression of the form:

$$
\begin{equation}
  \Pi(c, \eta, \nabla  \eta) = \int_{\Omega}    \left( f_{\alpha}(1-H) + f_{\beta}H \right)  + \frac{1}{2} | \gamma( \mathbf{n} ) \nabla  \eta |^2  + \frac{\delta^2}{2} (\Delta \eta)^2 ~dV
\end{equation}
$$

where $f_{\alpha}$ and $f_{\beta}$ are the free energy densities corresponding to $\alpha$ and $\beta$ phases, respectively, and are functions of composition $c$. $H$ is a function of the structural order parameter $\eta$.  $\delta$ is a scalar regularization parameter.  The interface normal vector $\mathbf{n}$ is given by

$$
\begin{equation}
\mathbf{n} = \frac{\nabla \eta}{|\nabla \eta|}
\end{equation}
$$

for $\nabla \eta \ne \mathbf{0}$, and $\mathbf{n} = \mathbf{0}$ when $\nabla \eta = \mathbf{0}$.

## Variational treatment
Following standard variational arguments (see Cahn-Hilliard formulation), we obtain the chemical potentials:

$$
\begin{align}
  \mu_{c}  &= (f_{\alpha,c}(1-H)+f_{\beta,c}H)
\end{align}
$$

$$
\begin{align}
  \mu_{\eta}  &= (f_{\beta,c}-f_{\alpha,c})H_{,\eta} - \nabla \cdot \mathbf{m} + \delta^2 \Delta(\Delta \eta)
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
  \frac{\partial c}{\partial t} &= -~\nabla \cdot (-M_c\nabla \mu_c)
\end{align}
$$

$$
\begin{align}
  &=M_c~\nabla \cdot (\nabla (f_{\alpha,c}(1-H)+f_{\beta,c}H))
  \end{align}
$$

and the PDE for Allen-Cahn dynamics is given by:

$$
  \begin{align}
    \frac{\partial \eta}{\partial t} &= -M_\eta \mu_\eta
\end{align}
$$

$$
\begin{align}
  &=-M_\eta ~ [(f_{\beta,c}-f_{\alpha,c})H_{,\eta} - \nabla \cdot \mathbf{m} + \delta^2 \Delta(\Delta \eta)]
\end{align}
$$

where $M_c$ and $M_\eta$ are the constant mobilities.  In order that the formulation only includes second order derivatives, an auxiliary field $\phi$ is introduced to break up the biharmonic term:

$$
\begin{align}
\phi = \Delta \eta
\end{align}
$$

and the PDE for Allen-Cahn dynamics becomes

$$
\begin{align}
    \frac{\partial \eta}{\partial t} =-M_\eta ~ ((f_{\beta,c}-f_{\alpha,c})H_{,\eta} - \nabla \cdot \mathbf{m}) + \delta^2 \Delta \phi .
\end{align}
$$

## Time discretization
Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:

$$
\begin{align}
 \phi^{n+1} &= \Delta \eta^n
\end{align}
$$

$$
\begin{align}
 \eta^{n+1} &= \eta^{n}  - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n -  \nabla \cdot \mathbf{m}^n + \delta^2 \Delta \phi^n)
\end{align}
$$

$$
\begin{align}
c^{n+1} &= c^{n}  + \Delta t M_{\eta}~\nabla \cdot (\nabla (f_{\alpha,c}^n(1-H^{n})+f_{\beta,c}^n H^{n}))
\end{align}
$$

## Weak formulation
In the weak formulation, considering an arbitrary variation $w$, the above equations can be expressed as residual equations.

$$
\begin{align}
  \int_{\Omega}   w  \phi^{n+1}  ~dV &=\int_{\Omega}  \nabla w \cdot \nabla \eta^n  ~dV
\end{align}
$$

$$
\begin{align}
r_{\phi x} &=  \nabla \eta^n
\end{align}
$$

$$
\begin{align}
  \int_{\Omega}   w  \eta^{n+1}  ~dV &= \int_{\Omega}   w \eta^{n} -   w    \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n - \kappa \Delta \eta^n)  ~dV
\end{align}
$$

$$
\begin{align}
  &=\int_{\Omega}  w  \left( \eta^{n} - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n) \right)+ \nabla w \cdot (- \Delta t M_{\eta})  ( \mathbf{m}^n - \delta^2 \phi^n) ~dV
\end{align}
$$

$$
\begin{align}
r_{\eta} &= \eta^{n} - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n)
\end{align}
$$

$$
\begin{align}
r_{\eta x} &= (- \Delta t M_{\eta})  ( \mathbf{m}^n - \delta^2 \phi^n)
\end{align}
$$

and

$$
\begin{align}
\int_{\Omega}   w  c^{n+1}  ~dV &= \int_{\Omega}   w c^{n} + w    \Delta t M_{c}~ \nabla \cdot (\nabla (f_{\alpha,c}^n(1-H^{n})+f_{\beta,c}^n H^{n}))  ~dV
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}   w c^{n} +  \nabla w   (-\Delta t M_{c})~ [~(f_{\alpha,cc}^n(1-H^{n})+f_{\beta,cc}^n H^{n}) \nabla c + ~((f_{\beta,c}^n-f_{\alpha,c}^n)H^{n}_{,\eta} \nabla \eta^n) ]  ~dV
\end{align}
$$

$$
\begin{align}
r_c &= c^{n}
\end{align}
$$

$$
\begin{align}
r_{cx} &= (-\Delta t M_{c})~ [~(f_{\alpha,cc}^n(1-H^{n})+f_{\beta,cc}^n H^{n}) \nabla c + ~((f_{\beta,c}^n-f_{\alpha,c}^n)H^{n}_{,\eta} \nabla \eta^n) ]
\end{align}
$$

The above values of $r_{\phi x}$, $r_{\eta}$, $r_{\eta x}$, $r_{c}$ and $r_{cx}$ are used to define the residuals in the following equations file: `textit{applications/anisotropyFacet/equations.h`

## Anisotropy
The above formulation is generic to any $\gamma(\mathbf{n})$.  In this application, we use an anisotropy of the form

$$
\begin{align}
\gamma(\mathbf{n}) = \gamma_0 \left( 1 - \sum_{i=1} \alpha_i (\mathbf m_i \cdot \mathbf n)^{w_i} \Theta(\mathbf m_i \cdot \mathbf n) \right),
\end{align}
$$

where $\mathbf{m}$ is a unit vector corresponding to a crystallographic orientation, $\gamma_0$ is a scaling factor for interfacial energy, $\alpha_i$ and $w_i$ are scalar parameters specific to each orientation, and $\Theta(\cdot)$ is the Heaviside function.  The derivatives with respect to components of the normal are

$$
\begin{align}
\frac{\partial \gamma(\mathbf{n})}{\partial n_j} = - \gamma_0  \sum_{i=1} w_i \alpha_i m_{ij} (\mathbf m_i \cdot \mathbf n)^{w_i-1} \Theta(\mathbf m_i \cdot \mathbf n) ,
\end{align}
$$

Calculation of $\gamma(\mathbf n)$ and  $\partial \gamma(\mathbf n)/\partial n_j$ is performed in an application-specific function located in `applications/anisotropyFacet/facet\_anisotropy.h`

This anisotropy was developed by M. Salvalaglio et al. (doi: 10.1021/acs.cgd.5b00165), and is extensively documented in their paper.  Briefly, we note that $\alpha_i$ determines the interfacial energy at the orientation $\mathbf m_i$, and $w_i$ determines how localized the change interfacial energy is around $\mathbf m_i$.  The Heaviside function $\Theta(\mathbf m_i \cdot \mathbf n)$ , which returns zero if $\mathbf m_i \cdot \mathbf n < 0$ and one otherwise, ensures that orientations are considered independently; i.e. there is no change in $\gamma(\mathbf n)$ around $-\mathbf m_1$ unless that corresponds to another listed orientation $\mathbf m_i$.  In its intended configuration, with $0 <\alpha_i < 1$ and high $w_i$ (e.g. $w_i = 50$), this anisotropy results in nearly flat facets at the orientations $\mathbf m_i$.
