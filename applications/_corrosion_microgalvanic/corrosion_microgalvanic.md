# Microgalvanic Corrosion (July 8, 2022)

This application simulates the evolution of the metal-electrolyte interface due to the microgalvanic coupling between the anodic and cathodic metals when immersed in the electrolyte. This model uses phase-field and smooth boundary methods to track the moving electrolyte/metal interface of the anodic phase. The electrolyte in this model is considered well-mixed; therefore, the model does not consider the effect of diffusion in the electrolyte. See Table 1 for the variable/parameter names in the code corresponding to each of the terms defined below.

## Free Energy Functional

Consider the free energy functional of the form

$$
\begin{equation}
 \widetilde {\mathcal F} = \int_\Omega \left(f_0\left(\{\varphi_j\},\psi\right)+\frac{\widetilde\epsilon^2}{2}\left(\sum^2_{j=1}|\nabla\varphi_j|^2+|\nabla\psi|^2\right)\right)d\Omega,
\end{equation}
$$

where $\widetilde {\mathcal F}$ represents the total free energy of the system scaled by the bulk free energy density coefficient, $W$, and $\tilde{\epsilon}^2=\epsilon/W$ is the rescaled gradient energy coefficient. Here  $\psi$ and $\varphi_j$  represent the order parameter for the electrolyte phase and for the $j^{th}$ metal phase respectively. The dimensionless bulk free energy density, $f_0$, is given by

$$
\begin{equation}
f_0\left(\{\varphi_j\},\psi\right)= \sum^N_{j=1}\left(\frac{\varphi_j^4}{4}-\frac{\varphi_j^2}{2}\right)+\left(\frac{\psi^4}{4}-\frac{\psi^2}{2}\right)+\frac{3}{2}\sum^N_{j=1}\sum^N_{k>j}\varphi_j^2\varphi_k^2+\frac{3}{2}\sum^N_{j=1} \varphi_j^2\psi^2.
\end{equation}
$$

In this application, we only consider two phases in the metal: the anodic phase, $\varphi_1$, and the cathodic phase, $\varphi_2$.


## Phase-field Equations
The evolution of the anodic metal/electrolyte interface is governed by Cahn-Hilliard equations with a source term for $\varphi_1$ and $\psi$,

$$
\begin{equation}
\frac{\partial \varphi_1}{\partial t} = \nabla \cdot \left(M\nabla\frac{\delta \widetilde F}{\delta \varphi_1} \right)+v\big|\nabla \psi \big|
\end{equation}
$$

and

$$
\begin{equation}
\frac{\partial \psi}{\partial t} = \nabla \cdot \left(M\nabla\frac{\delta \widetilde F}{\delta \psi} \right)-v\big|\nabla \psi \big|.
\end{equation}
$$

Since we do not consider any deposition on the cathodic phase, the evolution of $\varphi_2$ is given by the Cahn-Hilliard equation,

$$
\begin{equation}
\frac{\partial\varphi_2}{\partial t}=\nabla \cdot \left(M\nabla\frac{\delta \widetilde F}{\delta \varphi_2} \right)
\end{equation}
$$

In the previous Eqs., $M$ is the Cahn-Hilliard mobility coefficient, assumed to be equal for all phases, and $v$ is the velocity of the interface normal to the surface, which is related to the reaction current at the interface via Faraday’s law of electrolysis:

$$
\begin{equation}
v = -\frac{V_m \xi_1 i_{rxn,1}}{z_mF}.
\end{equation}
$$

Here $V_m$ is the molar volume of the metal, $i_{rxn,1}$, is the anodic current density, $z_m$ is the dissolved metal cation charge number, and $F$ is Faraday’s constant. The variable $\xi_1$ is used as a weighing factor for the contribution of the anodic current density over the total current density. This factor is defined in the following section. The interfacial mobility, $M$, is set as a function of the local current density,

$$
\begin{equation}
M = 2\delta\frac{V_m|\xi_1 i_{rxn,1}|\psi}{z_mF},
\end{equation}
$$

where $2\delta=2\sqrt{2\widetilde{\epsilon}^2}$ gives the equilibrium interfacial thickness.

## Electrochemical Equations

In general, the transport of the ionic species in the electrolyte is controlled by migration, diffusion, and convection. Though PRISMS-PF is capable of including diffusion, in this application, the electrolyte is considered to be well-stirred. The current density ${i_e}$ in the electrolyte is given by

$$
\begin{equation}
{i_e} = -\kappa_e\nabla\Phi_e
\end{equation}
$$

where $\kappa_e$ is the conductivity of the electrolyte while $\Phi_e$ is the electrostatic potential in the electrolyte. To satisfy the current continuity in the electrolyte,

$$
\begin{equation}
\nabla \cdot {i_e} = 0.
\end{equation}
$$

At the metal/electrolyte interface, the normal component of ${i_e}$ is equal to $i_{rxn}$

$$
\begin{equation}
{i_e}\cdot{\hat n_{m/e}}=i_{rxn}.
\end{equation}
$$

where ${\hat n_{m/e}}$ is the unit normal vector of the metal/electrolyte interface pointing out of the metal. The previous Eqs. can combined using the Smoothed Boundary Method (SBM) [H.-C. Yu, H.-Y. Chen, and K. Thornton, Model. Simul. Mater. Sci. Eng. **20**, 075008
(2012)] to obtain an expression for $\Phi_e$ within the electrolyte that includes the appropriate boundary condition at the metal/electrolyte interface. For this system, the SBM is applied with $\psi$ as the domain parameter. The resulting equation for $\Phi_e$ is given by

$$
\begin{equation}
\nabla\cdot\left(\psi\kappa_e\nabla\Phi_e\right) =-|\nabla\psi|i_{rxn}.
\end{equation}
$$

The anodic current density is given by

$$
\begin{equation}
i_{rxn,1} = \frac{i_{corr,1}\cdot e^{\frac{\eta_1}{A_1}}}{1+\frac{i_{corr,1}\cdot e^{\frac{\eta_1}{A_{1}}}}{i_{max}}},
\end{equation}
$$

where $i_{corr,1}$ is the corrosion current density of the anodic phase, $\eta_1$ is the anodic overpotential, given by

$$
\begin{equation}
\eta_1 = \Phi_m - \Phi_e - E_{corr,1},
\end{equation}
$$

where $\Phi_m$ and $E_{corr,1}$ represent the applied electrostatic potential (we set is as zero because we consider free immersion) and the corrosion potential of the anodic metal, respectively. The constant $A_1$ is the Tafel slope of the anodic metal and $i_{max}$ is the maximum current density. The cathodic current density is given as

$$
\begin{equation}
i_{rxn,2}=-i_{corr,2}\cdot e^{\frac{\eta_2}{A_2}},
\end{equation}
$$

where $i_{corr,2}$ is the corrosion current density of the cathodic phase, $\eta_2$ is the cathodic overpotential which is given in the same form as its anodic counterpart:

$$
\begin{equation}
\eta_2 = \Phi_m - \Phi_e - E_{corr,2},
\end{equation}
$$

where $E_{corr,2}$ represents the corrosion potential of the cathodic metal. The reaction current density $i_{rxn}$ is formulated as a linear combination of $i_{rxn,1}$ and $i_{rxn,2}$

$$
\begin{equation}
i_{rxn} = \sum_{j=1}^{2}\xi_{j}i_{rxn,j},
\end{equation}
$$

where $\xi_j$ is an interpolation function that varies from 0 to 1, which indicates the region of each of the metal phases. We define $\xi$ for each of the metal phases as

$$
\begin{equation}
\xi_{2} = \frac{\varphi_{2}}{max\left(\varphi_1 + \varphi_2,\zeta\right)},
\end{equation}
$$

and

$$
\begin{equation}
\xi_{1} = 1- \varphi_2,
\end{equation}
$$

where $\zeta$ is small regularization constant  added to the denominator to avoid division by zero.

## Time Discretization

The fields $\varphi_1$, $\varphi_2$ and $\psi$ are solved using an explicit Euler method for time integration. For we employ a splitting strategy to transform Eqs. 3, 4, and 5, into second order PDEs

$$
\begin{equation}
\frac{\partial \varphi_1}{\partial t} = \nabla \cdot \left(M\nabla\mu_1 \right)+v\big|\nabla \psi \big|,
\end{equation}
$$

$$
\begin{equation}
\frac{\partial\varphi_2}{\partial t}=\nabla \cdot \left(M\nabla\mu_2 \right),
\end{equation}
$$

and

$$
\begin{equation}
\frac{\partial \psi}{\partial t} = \nabla \cdot \left(M\nabla\mu_{\psi} \right)-v\big|\nabla \psi \big|,
\end{equation}
$$

where

$$
\begin{equation}
\mu_{1}=f_1-\widetilde{\epsilon}^2\nabla^2\varphi_{1},
\end{equation}
$$

$$
\begin{equation}
\mu_{2}=f_2-\widetilde{\epsilon}^2\nabla^2\varphi_{2},
\end{equation}
$$

and

$$
\begin{equation}
\mu_{\psi}=f_{\psi}-\widetilde{\epsilon}^2\nabla^2\psi.
\end{equation}
$$

The terms $f_1$, $f_2$ and $f_\psi$ are given by

$$
\begin{equation}
f_1=W\left(\varphi_{1}^3-\varphi_{1}+3\varphi_{1}\varphi_{2}^2+3\varphi_{1}\psi^2\right),
\end{equation}
$$

$$
\begin{equation}
f_2=W\left(\varphi_{2}^3-\varphi_{2}+3\varphi_{2}\varphi_{1}^2+3\varphi_{2}\psi^2\right).
\end{equation}
$$

and

$$
\begin{equation}
f_{\psi}= W\left(\psi^3-\psi+3\psi\varphi_{1}^2+3\psi\varphi_{2}^2\right),
\end{equation}
$$

Considering forward Euler explicit time stepping, the time-discretized versions of the equations above are

$$
\begin{equation}
\varphi^{n+1}_1 = \varphi^n_1+\Delta t \left(\nabla\cdot \left(M\nabla\mu^n_1\right)+v|\nabla \psi^n|\right),
\end{equation}
$$

$$
\begin{equation}
\varphi^{n+1}_2 = \varphi^n_2+\Delta t \nabla\cdot \left(M\nabla\mu^n_2\right),
\end{equation}
$$

and

$$
\begin{equation}
\psi^{n+1} = \psi^n+\Delta t \left(\nabla\cdot\left(M\nabla\mu^n_{\psi}\right)-v|\nabla\psi^n|\right),
\end{equation}
$$

where

$$
\begin{equation}
\mu^n_1=f^n_1-\widetilde{\epsilon}^2\nabla^2\varphi^n_1,
\end{equation}
$$

$$
\begin{equation}
\mu^n_2=f^n_2-\widetilde{\epsilon}^2\nabla^2\varphi^n_2,
\end{equation}
$$

and

$$
\begin{equation}
\mu^n_{\psi}=f^n_{\psi}-\widetilde{\epsilon}^2\nabla^2\psi^n.
\end{equation}
$$

## Weak Formulation

The weak formulation for Eqs. 28-33 is derived by multiplying by a test function, $\omega$, and integrating over the whole system. After applying the divergence theorem and rearranging in order to obtain only first-order spatial derivatives, the terms denoted as 'value' terms are those that multiply the test function and the terms denoted as 'gradient' terms are those that multiply by the gradient of the test function, $\nabla \omega$(see `https://prisms-center.github.io/phaseField/doxygen_files/app_files.html` for details). The 'value' and 'gradient' terms are obtained from the following equations:

$$
\begin{align}
\int_\Omega \omega \varphi_1^{n+1}dV = \int_\Omega \omega\left(\varphi_1^{n}+\Delta tv|\nabla\psi^n|\right)dV + \int_\Omega \nabla\omega\cdot \left(-\Delta t M\nabla\mu_1^{n}\right) dV
\end{align}
$$

$$
\begin{align}
r_1 &= \left(\varphi^{n}_1+\Delta tv|\nabla\psi^n|\right)
\end{align}
$$

$$
\begin{align}
r_{1 x} &= \left(-\Delta t M\nabla\mu_1^{n}\right)
\end{align}
$$

$$
\begin{align}
\int_\Omega \omega \varphi_2^{n+1}dV = \int_\Omega \omega \varphi_2^n dV + \int_\Omega \nabla\omega\cdot \left(-\Delta t M\nabla\mu_2^{n}\right) dV
\end{align}
$$

$$
\begin{align}
r_2 &= \varphi^n_2
\end{align}
$$

$$
\begin{align}
r_{2x} &= \left(-\Delta t M\nabla\mu^{n}_2\right)
\end{align}
$$

$$
\begin{align}
\int_\Omega \omega \psi^{n+1}dV = \int_\Omega \omega \left(\psi^n-\Delta tv|\nabla\psi^n|\right) dV + \int_\Omega \nabla\omega\cdot \left(-\Delta t M\nabla\mu^n_{\psi}\right) dV
\end{align}
$$

$$
\begin{align}
r_{\psi} &= \left(\psi^n-\Delta tv|\nabla\psi^n|\right)
\end{align}
$$

$$
\begin{align}
r_{\psi x} &= \left(-\Delta t M\nabla\mu^n_{\psi}\right)
\end{align}
$$

$$
\begin{align}
\int_\Omega \omega\mu_1^{n+1}dV = \int_\Omega \omega f^n_1 dV +\int_\Omega \nabla\omega\cdot \left(\widetilde{\epsilon}^2\nabla\varphi^n_1\right) dV
\end{align}
$$

$$
\begin{align}
r_{\mu 1} &= f^n_1
\end{align}
$$

$$
\begin{align}
r_{\mu 1 x} &= \left(\widetilde{\epsilon}^2\nabla\varphi^n_1\right)
\end{align}
$$

$$
\begin{align}
\int_\Omega \omega\mu_2^{n+1} dV = \int_\Omega \omega f^n_2 dV +\int_\Omega \nabla\omega\cdot \left(\widetilde{\epsilon}^2\nabla\varphi^n_2\right) dV
\end{align}
$$

$$
\begin{align}
r_{\mu 2} &= f^n_2
\end{align}
$$

$$
\begin{align}
r_{\mu 2 x} &= \left(\widetilde{\epsilon}^2\nabla\varphi^n_2\right)
\end{align}
$$

$$
\begin{align}
\int_\Omega \omega\mu_{\psi}^{n+1}dV = \int_\Omega \omega f^n_{\psi} dV +\int_\Omega \nabla\omega\cdot \left(\widetilde{\epsilon}^2\nabla\psi^n\right) dV
\end{align}
$$

$$
\begin{align}
r_{\mu\psi} &= f^n_{\psi}
\end{align}
$$

$$
\begin{align}
r_{\mu\psi x} &= \left(\widetilde{\epsilon}^2\nabla\psi^n\right)
\end{align}
$$

For the weak formulation of the time-independent electrostatic potential equation we need to specify LHS and RHS terms. To write the equation in terms of a Newton iteration, the solution, $\Phi_e$, can be written as the sum of an initial guess, $\Phi^0_e$, and an update, $\Delta \Phi_e$ (see `https://prisms-center.github.io/phaseField/doxygen_files/app_files.html` for details):

$$
\begin{align}
\int_\Omega \omega \left(\left(\frac{\partial i_{rxn}}{\partial\Phi_e}\right)_{\Phi_e = \Phi_e^0}|\nabla\psi|\Delta\Phi_e\right)dV
\end{align}
$$

$$
\begin{align}
&+ \int_\Omega \nabla\omega \cdot \left(-\psi\kappa_e\nabla\left(\Delta\Phi_e\right)\right)dV=\int_\Omega\omega \left(-|\nabla\psi|i_{rxn}\right)dV+\int_\Omega\nabla\omega\cdot \left(\psi\kappa_e\nabla\Phi^0_e\right)dV,
\end{align}
$$

$$
\begin{align}
r_{\Delta \Phi_e} &= \left(\left(\frac{\partial i_{rxn}}{\partial\Phi_e}\right)_{\Phi_e = \Phi^0_e}|\nabla\psi|\Delta\Phi_e\right)
\end{align}
$$

$$
\begin{align}
r_{\Delta \Phi_e x} &= \left(-\psi\kappa_e\nabla\left(\Delta\Phi_e\right)\right)
\end{align}
$$

$$
\begin{align}
r_{\Phi^0_e} &= \left(-|\nabla\psi|i_{rxn}\right)
\end{align}
$$

$$
\begin{align}
r_{\Phi^0_ex} &= \left(\psi\kappa_e\nabla\Phi^0_e\right)
\end{align}
$$

where

$$
\begin{align}
\frac{\partial i_{rxn}}{\partial \Phi_e} = -\xi_{1}\left(\frac{i_{max}}{i_{max}+i_{corr,1}\cdot e^{\frac{\eta_1}{A_{1}}}}\right)^2\left(\frac{i_{corr,1}\cdot e^{\frac{\eta_1}{A_{1}}}}{A_1}\right)-\xi_{2}\left(\frac{i_{rxn,2}}{A_2}\right).
\end{align}
$$

The above values of are used to define the equation terms in the input file:
`applications/corrosion\_microgalvanic/equations.cc`





 |Name in the model equations | Variable/parameter name in the code |
 |----------------------------|-------------------------------------|
| $\varphi_1$           |     nAnodic              |
| $\varphi_2$           |   nCathodic              |
| $\psi$                |    psi                   |
| $v$                    |    v                    |
| $M$                    |    MnV                  |
| $V_m$                  |    VMV                  |
| $z_m$                  |    zMV                  |
| $\delta$                |   delta                |
| $\widetilde{{\epsilon}}^2$ | epssqV              |
| $\kappa_e$            |     kappa                |
| $\Phi_e$              |     Phi                  |
| $i_{corr,1}$          |     i0Anodic             |
| $i_{corr,2}$          |     i0Cathodic           |
| $i_{rxn,1}$           |     iAnodic              |
| $i_{rxn,2}$           |     iCathodic            |
| $i_{max}$             |     iMax                 |
| $\eta_1$              |     etaAnodic            |
| $\eta_2$              |     etaCathodic          |
| \xi_1$                |    xiAnodic              |
| $\xi_2$               |     xiCathodic           |
| $A_1$                 |     AAnodic              |
| $A_2$                 |     ACathodic            |
| $E_{corr,1}$          |     EcorrAnodic          |
| $E_{corr,2}$          |     EcorrCathodic        |
| $\zeta$               |     lthresh              |


Table 1: Variables/parameters names used in the model equations and the corresponding names used in the code
