# Microgalvanic Corrosion (July 8, 2022)

This application simulates the evolution of the metal-electrolyte interface due to the microgalvanic coupling between the anodic and cathodic metals when immersed in the electrolyte. This model uses phase-field and smooth boundary methods to track the moving electrolyte/metal interface of the anodic phase. The electrolyte in this model is considered well-mixed; therefore, the model does not consider the effect of diffusion in the electrolyte. See Table 1 for the variable/parameter names in the code corresponding to each of the terms defined below.

\section{Free Energy Functional}

Consider the free energy functional of the form

\begin{equation} \label{free_energy}
 \widetilde {\mathcal F} = \int_\Omega \left(f_0\left(\{\varphi_j\},\psi\right)+\frac{\widetilde\epsilon^2}{2}\left(\sum^2_{j=1}|\nabla\varphi_j|^2+|\nabla\psi|^2\right)\right)d\Omega,
\end{equation}

where $\widetilde {\mathcal F}$ represents the total free energy of the system scaled by the bulk free energy density coefficient, $W$, and $\tilde{\epsilon}^2=\epsilon/W$ is the rescaled gradient energy coefficient. Here  $\psi$ and $\varphi_j$  represent the order parameter for the electrolyte phase and for the $j^{th}$ metal phase respectively. The dimensionless bulk free energy density, $f_0$, is given by

\begin{equation} \label{f0}
f_0\left(\{\varphi_j\},\psi\right)= \sum^N_{j=1}\left(\frac{\varphi_j^4}{4}-\frac{\varphi_j^2}{2}\right)+\left(\frac{\psi^4}{4}-\frac{\psi^2}{2}\right)+\frac{3}{2}\sum^N_{j=1}\sum^N_{k>j}\varphi_j^2\varphi_k^2+\frac{3}{2}\sum^N_{j=1} \varphi_j^2\psi^2.
\end{equation}

In this application, we only consider two phases in the metal: the anodic phase, $\varphi_1$, and the cathodic phase, $\varphi_2$.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bigskip
\section{Phase-field Equations}
The evolution of the anodic metal/electrolyte interface is governed by Cahn-Hilliard equations with a source term for $\varphi_1$ and $\psi$,

\begin{equation} \label{dphi1dt}
\frac{\partial \varphi_1}{\partial t} = \nabla \cdot \left(M\nabla\frac{\delta \widetilde \mathcal F}{\delta \varphi_1} \right)+v\big|\nabla \psi \big|
\end{equation}

and

\begin{equation} \label{dpsidt}
\frac{\partial \psi}{\partial t} = \nabla \cdot \left(M\nabla\frac{\delta \widetilde \mathcal F}{\delta \psi} \right)-v\big|\nabla \psi \big|.
\end{equation}

Since we do not consider any deposition on the cathodic phase, the evolution of $\varphi_2$ is given by the Cahn-Hilliard equation, 

\begin{equation} \label{dphi2dt}
\frac{\partial\varphi_2}{\partial t}=\nabla \cdot \left(M\nabla\frac{\delta \widetilde \mathcal F}{\delta \varphi_2} \right)
\end{equation}

In Eqs. \eqref{dphi1dt}, \eqref{dpsidt}, and \eqref{dphi2dt}, $M$ is the Cahn-Hilliard mobility coefficient, assumed to be equal for all phases, and $v$ is the velocity of the interface normal to the surface, which is related to the reaction current at the interface via Faraday’s law of electrolysis:

\begin{equation} \label{int_vel}
v = -\frac{V_m \xi_1 i_{rxn,1}}{z_mF}.
\end{equation}

Here $V_m$ is the molar volume of the metal, $i_{rxn,1}$, is the anodic current density, $z_m$ is the dissolved metal cation charge number, and $F$ is Faraday’s constant. The variable $\xi_1$ is used as a weighing factor for the contribution of the anodic current density over the total current density. This factor is defined in the following section. The interfacial mobility, $M$, is set as a function of the local current density,

\begin{equation} \label{mobility}
M = 2\delta\frac{V_m|\xi_1 i_{rxn,1}|\psi}{z_mF},
\end{equation}

where $2\delta=2\sqrt{2\widetilde{\epsilon}^2}$ gives the equilibrium interfacial thickness.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bigskip
\section{Electrochemical Equations}

In general, the transport of the ionic species in the electrolyte is controlled by migration, diffusion, and convection. Though PRISMS-PF is capable of including diffusion, in this application, the electrolyte is considered to be well-stirred. The current density $\vb*{i_e}$ in the electrolyte is given by

\begin{equation} \label{electrolyte_i_e}
\vb*{i_e} = -\kappa_e\nabla\Phi_e
\end{equation}

where $\kappa_e$ is the conductivity of the electrolyte while $\Phi_e$ is the electrostatic potential in the electrolyte. To satisfy the current continuity in the electrolyte,

\begin{equation} \label{curr_continuity}
\nabla \cdot \vb*{i_e} = 0.
\end{equation}

At the metal/electrolyte interface, the normal component of $\vb*{i_e}$ is equal to $i_{rxn}$

\begin{equation} \label{i_e_BC}
\vb*{i_e}\cdot\vb*{\hat n_{m/e}}=i_{rxn}.
\end{equation}

where $\vb*{\hat n_{m/e}}$ is the unit normal vector of the metal/electrolyte interface pointing out of the metal. \Cref{electrolyte_i_e,curr_continuity,i_e_BC} can combined using the Smoothed Boundary Method (SBM) [H.-C. Yu, H.-Y. Chen, and K. Thornton, Model. Simul. Mater. Sci. Eng. {\bf 20}, 075008
(2012)] to obtain an expression for $\Phi_e$ within the electrolyte that includes the appropriate boundary condition at the metal/electrolyte interface. For this system, the SBM is applied with $\psi$ as the domain parameter. The resulting equation for $\Phi_e$ is given by 

\begin{equation} \label{Phi_e}
\nabla\cdot\left(\psi\kappa_e\nabla\Phi_e\right) =-|\nabla\psi|i_{rxn}.
\end{equation}

The anodic current density is given by

\begin{equation} \label{i_rxn_1}
i_{rxn,1} = \frac{i_{corr,1}\cdot e^{\frac{\eta_1}{A_1}}}{1+\frac{i_{corr,1}\cdot e^{\frac{\eta_1}{A_{1}}}}{i_{max}}},
\end{equation}

where $i_{corr,1}$ is the corrosion current density of the anodic phase, $\eta_1$ is the anodic overpotential, given by

\begin{equation} \label{eta_1}
\eta_1 = \Phi_m - \Phi_e - E_{corr,1},
\end{equation}

where $\Phi_m$ and $E_{corr,1}$ represent the applied electrostatic potential (we set is as zero because we consider free immersion) and the corrosion potential of the anodic metal, respectively. The constant $A_1$ is the Tafel slope of the anodic metal and $i_{max}$ is the maximum current density. The cathodic current density is given as

\begin{equation} \label{i_rxn_2}
i_{rxn,2}=-i_{corr,2}\cdot e^{\frac{\eta_2}{A_2}},
\end{equation}

where $i_{corr,2}$ is the corrosion current density of the cathodic phase, $\eta_2$ is the cathodic overpotential which is given in the same form as its anodic counterpart:

\begin{equation} \label{eta_2}
\eta_2 = \Phi_m - \Phi_e - E_{corr,2},
\end{equation}

where $E_{corr,2}$ represents the corrosion potential of the cathodic metal. The reaction current density $i_{rxn}$ is formulated as a linear combination of $i_{rxn,1}$ and $i_{rxn,2}$

\begin{equation} \label{i_rxn}
i_{rxn} = \sum^{2}_{j=1}\xi_{j}i_{rxn,j},
\end{equation}

where $\xi_j$ is an interpolation function that varies from 0 to 1, which indicates the region of each of the metal phases. We define $\xi$ for each of the metal phases as

\begin{equation} \label{xi_2}
\xi_{2} = \frac{\varphi_{2}}{max\left(\varphi_1 + \varphi_2,\zeta\right)},
\end{equation}

and

\begin{equation} \label{xi_1}
\xi_{1} = 1- \varphi_2,
\end{equation}

where $\zeta$ is small regularization constant  added to the denominator to avoid division by zero.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bigskip
\section{Time Discretization}

The fields $\varphi_1$, $\varphi_2$ and $\psi$ are solved using an explicit Euler method for time integration. For we employ a splitting strategy to transform Eqs. \eqref{dphi1dt}, \eqref{dpsidt}, and \eqref{dphi2dt}, into second order PDEs

\begin{equation} \label{dphi1dt_split}
\frac{\partial \varphi_1}{\partial t} = \nabla \cdot \left(M\nabla\mu_1 \right)+v\big|\nabla \psi \big|,
\end{equation}

\begin{equation} \label{dphi2dt_split}
\frac{\partial\varphi_2}{\partial t}=\nabla \cdot \left(M\nabla\mu_2 \right),
\end{equation}
and
\begin{equation} \label{dpsidt_split}
\frac{\partial \psi}{\partial t} = \nabla \cdot \left(M\nabla\mu_{\psi} \right)-v\big|\nabla \psi \big|,
\end{equation}

where 

\begin{equation} \label{mu_phi1}
\mu_{1}=f_1-\widetilde{\epsilon}^2\nabla^2\varphi_{1},
\end{equation}

\begin{equation} \label{mu_phi2}
\mu_{2}=f_2-\widetilde{\epsilon}^2\nabla^2\varphi_{2},
\end{equation}
and
\begin{equation} \label{mu_psi}
\mu_{\psi}=f_{\psi}-\widetilde{\epsilon}^2\nabla^2\psi.
\end{equation}

The terms $f_1$, $f_2$ and $f_\psi$ are given by

\begin{equation} \label{f_1}
f_1=W\left(\varphi_{1}^3-\varphi_{1}+3\varphi_{1}\varphi_{2}^2+3\varphi_{1}\psi^2\right),
\end{equation}

\begin{equation} \label{f_2}
f_2=W\left(\varphi_{2}^3-\varphi_{2}+3\varphi_{2}\varphi_{1}^2+3\varphi_{2}\psi^2\right).
\end{equation}
and
\begin{equation} \label{f_psi}
f_{\psi}= W\left(\psi^3-\psi+3\psi\varphi_{1}^2+3\psi\varphi_{2}^2\right),
\end{equation}

Considering forward Euler explicit time stepping, the time-discretized versions of the equations above are

\begin{equation} \label{phi1_np1}
\varphi^{n+1}_1 = \varphi^n_1+\Delta t \left(\nabla\cdot \left(M\nabla\mu^n_1\right)+v|\nabla \psi^n|\right),
\end{equation}

\begin{equation} \label{phi2_np1}
\varphi^{n+1}_2 = \varphi^n_2+\Delta t \nabla\cdot \left(M\nabla\mu^n_2\right),
\end{equation}
and
\begin{equation} \label{psi_np1}
\psi^{n+1} = \psi^n+\Delta t \left(\nabla\cdot\left(M\nabla\mu^n_{\psi}\right)-v|\nabla\psi^n|\right),
\end{equation}
where
\begin{equation} \label{mu1_n}
\mu^n_1=f^n_1-\widetilde{\epsilon}^2\nabla^2\varphi^n_1,
\end{equation}

\begin{equation} \label{mu2_n}
\mu^n_2=f^n_2-\widetilde{\epsilon}^2\nabla^2\varphi^n_2,
\end{equation}
and
\begin{equation} \label{mu_psi_n}
\mu^n_{\psi}=f^n_{\psi}-\widetilde{\epsilon}^2\nabla^2\psi^n.
\end{equation}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bigskip
\section{Weak Formulation}

The weak formulation for Eqs. \labelcref{phi1_np1,phi2_np1,psi_np1,mu1_n,mu2_n,mu_psi_n} is derived by multiplying by a test function, $\omega$, and integrating over the whole system. After applying the divergence theorem and rearranging in order to obtain only first-order spatial derivatives, the terms denoted as 'value' terms are those that multiply the test function and the terms denoted as 'gradient' terms are those that multiply by the gradient of the test function, $\nabla \omega$(see \href{https://prisms-center.github.io/phaseField/doxygen_files/app_files.html}{User Manual - equations.cc} for details). The 'value' and 'gradient' terms are obtained from the following equations:

\begin{equation} \label{phi1_weak}
\int_\Omega \omega \varphi^{n+1}_1dV = \int_\Omega \omega\underbrace{\left(\varphi^{n}_1+\Delta tv|\nabla\psi^n|\right)}_{r_1}dV + \int_\Omega \nabla\omega\cdot\underbrace{\left(-\Delta t M\nabla\mu^{n}_1\right)}_{r_{1 x}}dV
\end{equation}

\begin{equation} \label{phi2_weak}
\int_\Omega \omega \varphi^{n+1}_2dV = \int_\Omega \omega\underbrace{\varphi^n_2}_{r_2}dV + \int_\Omega \nabla\omega\cdot\underbrace{\left(-\Delta t M\nabla\mu^{n}_2\right)}_{r_{2x}}dV
\end{equation}

\begin{equation} \label{psi_weak}
\int_\Omega \omega \psi^{n+1}dV = \int_\Omega \omega\underbrace{\left(\psi^n-\Delta tv|\nabla\psi^n|\right)}_{r_{\psi}}dV + \int_\Omega \nabla\omega\cdot\underbrace{\left(-\Delta t M\nabla\mu^n_{\psi}\right)}_{r_{\psi x}}dV
\end{equation}

\begin{equation} \label{mu1_weak}
\int_\Omega \omega\mu^{n+1}_1dV = \int_\Omega \omega \underbrace{f^n_1}_{r_{\mu 1}} dV +\int_\Omega \nabla\omega\cdot\underbrace{\left(\widetilde{\epsilon}^2\nabla\varphi^n_1\right)}_{r_{\mu 1 x}}dV
\end{equation}

\begin{equation} \label{mu2_weak}
\int_\Omega \omega\mu^{n+1}_2dV = \int_\Omega \omega \underbrace{f^n_2}_{r_{\mu 2}} dV +\int_\Omega \nabla\omega\cdot\underbrace{\left(\widetilde{\epsilon}^2\nabla\varphi^n_2\right)}_{r_{\mu 2 x}}dV
\end{equation}

\begin{equation} \label{mu_psi_weak}
\int_\Omega \omega\mu^{n+1}_{\psi}dV = \int_\Omega \omega \underbrace{f^n_{\psi}}_{r_{\mu\psi}} dV +\int_\Omega \nabla\omega\cdot\underbrace{\left(\widetilde{\epsilon}^2\nabla\psi^n\right)}_{r_{\mu\psi x}}dV
\end{equation}

For the weak formulation of the time-independent electrostatic potential equation \labelcref{Phi_e} we need to specify LHS and RHS terms. To write the equation in terms of a Newton iteration, the solution, $\Phi_e$, can be written as the sum of an initial guess, $\Phi^0_e$, and an update, $\Delta \Phi_e$ (see \href{https://prisms-center.github.io/phaseField/doxygen_files/app_files.html}{User Manual - nonExplicitEquationRHS} for details):

\begin{equation} \label{Phi_weak}
\int _\Omega \omega\underbrace{\left(\left(\frac{\partial i_{rxn}}{\partial\Phi_e}\right)_{\Phi_e = \Phi^0_e}|\nabla\psi|\Delta\Phi_e\right)}_{r_{\Delta \Phi_e}}dV+\int_\Omega \nabla\omega \cdot\underbrace{\left(-\psi\kappa_e\nabla\left(\Delta\Phi_e\right)\right)}_{r_{\Delta \Phi_e x}}dV=\int_\Omega\omega\underbrace{\left(-|\nabla\psi|i_{rxn}\right)}_{r_{\Phi^0_e}}dV+\int_\Omega\nabla\omega\cdot\underbrace{\left(\psi\kappa_e\nabla\Phi^0_e\right)}_{r_{\Phi^0_ex}}dV,
\end{equation}

where

\begin{equation} \label{dirxndPhi}
\frac{\partial i_{rxn}}{\partial \Phi_e} = -\xi_{1}\left(\frac{i_{max}}{i_{max}+i_{corr,1}\cdot e^{\frac{\eta_1}{A_{1}}}}\right)^2\left(\frac{i_{corr,1}\cdot e^{\frac{\eta_1}{A_{1}}}}{A_1}\right)-\xi_{2}\left(\frac{i_{rxn,2}}{A_2}\right).
\end{equation}

The above values of are used to define the equation terms in the input file:\\
{\it applications/corrosion\_microgalvanic/equations.cc}
