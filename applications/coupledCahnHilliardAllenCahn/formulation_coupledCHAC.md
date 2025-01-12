# PRISMS PhaseField: Coupled Cahn-Hilliard and Allen-Cahn Dynamics

Consider a free energy expression of the form:
\begin{equation}
  \Pi(c, \eta, \grad  \eta) = \int_{\Omega}    \left( f_{\alpha}(1-H) + f_{\beta}H \right)  + \frac{\kappa}{2} \grad  \eta  \cdot \grad  \eta    ~dV 
\end{equation}
where $f_{\alpha}$ and $f_{\beta}$ are the free energy densities corresponding to $\alpha$ and $\beta$ phases, respectively, and are functions of composition $c$. $H$ is a function of the structural order parameter $\eta$. Note that we don't have the gradient terms for the composition, i.e, $\grad c$ terms, unlike classical Cahn-Hillard formulation.

\section{Variational treatment}
Following standard variational arguments (see Cahn-Hilliard formulation), we obtain the chemical potentials:
\begin{align}
  \mu_{c}  &= (f_{\alpha,c}(1-H)+f_{\beta,c}H)  \\
  \mu_{\eta}  &= (f_{\beta}-f_{\alpha})H_{,\eta} - \kappa \Delta \eta 
\end{align}

\section{Kinetics}
Now the PDE for Cahn-Hilliard dynamics is given by:
\begin{align}
  \frac{\partial c}{\partial t} &= -~\grad \cdot (-M_c\grad \mu_c)\\
  &=M_c~\grad \cdot (\grad (f_{\alpha,c}(1-H)+f_{\beta,c}H)) 
  \end{align}
  and the PDE for Allen-Cahn dynamics is given by:
  \begin{align}
    \frac{\partial \eta}{\partial t} &= -M_\eta \mu_\eta \\
  &=-M_\eta ~ ((f_{\beta}-f_{\alpha})H_{,\eta} - \kappa \Delta \eta) 
\end{align}
where $M_c$ and $M_\eta$ are the constant mobilities. 

\section{Time discretization}
Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:
\begin{align}
 \eta^{n+1} &= \eta^{n}  - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n - \kappa \Delta \eta^n) \\
c^{n+1} &= c^{n}  + \Delta t M_{\eta}~\grad \cdot (\grad (f_{\alpha,c}^n(1-H^{n})+f_{\beta,c}^n H^{n}))
\end{align}

\section{Weak formulation}
In the weak formulation, considering an arbitrary variation $w$, the above equations can be expressed as residual equations:
\begin{align}
  \int_{\Omega}   w  \eta^{n+1}  ~dV &= \int_{\Omega}   w \eta^{n} -   w    \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n - \kappa \Delta \eta^n)  ~dV\\
  &=\int_{\Omega}  w  \left( \underbrace{\eta^{n} - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n)}_{r_{\eta}} \right)+ \grad w \cdot \underbrace{(- \Delta t M_{\eta}\kappa) \grad \eta^{n}}_{r_{\eta x}} ~dV 
\end{align}
and 
\begin{align}
  \int_{\Omega}   w  c^{n+1}  ~dV &= \int_{\Omega}   w c^{n} + w    \Delta t M_{c}~ \grad \cdot (\grad (f_{\alpha,c}^n(1-H^{n})+f_{\beta,c}^n H^{n}))  ~dV\\
    &= \int_{\Omega}   w \underbrace{c^{n}}_{r_c} +  \grad w   \underbrace{(-\Delta t M_{c})~ [~(f_{\alpha,cc}^n(1-H^{n})+f_{\beta,cc}^n H^{n}) \grad c + ~((f_{\beta,c}^n-f_{\alpha,c}^n)H^{n}_{,\eta} \grad \eta) ] }_{r_{cx}} ~dV
\end{align}

\vskip 0.25in
The above values of $r_{\eta}$, $r_{\eta x}$, $r_{c}$ and $r_{cx}$ are used to define the residuals in the following parameters file: \\
\textit{applications/coupledCahnHilliardAllenCahn/parameters.h}
