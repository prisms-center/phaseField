# Corrosion (March 3, 2024)

This application simulates the evolution of the metal-electrolyte interface during the anodic corrosion reaction. The model employed  [Chadwick et al., J. Electrochem. Soc.,{\bf 10}, C633-C646 (2018)]  uses the phase-field and smoothed boundary methods to track the moving metal/electrolyte interface and to couple it to mass transport (diffusion and migration) within the electrolyte and Butler-Volmer electrochemical kinetics.

\section{Free Energy}
Consider a  free energy expression of the form:
\begin{equation}
  \Pi(\eta,\psi) = \int_{\Omega}  \left( W f_0(\eta,\psi)+ \frac{\epsilon^2}{2} |\grad \eta|^2 +  \frac{\epsilon^2}{2} |\grad \psi|^2 \right) dV
\end{equation}
\noindent with $f_0$ given by
\begin{equation}
 f_0(\eta,\psi) = \left( \frac{\eta^4}{4} - \frac{\eta^2}{2} \right) + \left( \frac{\psi^4}{4} - \frac{\psi^2}{2} \right) + \gamma ~\eta^2 \psi^2,
\end{equation}
\noindent where the fields $\eta$  and $\psi$ represent the physical domain of the metal and the liquid electrolyte, respectively, $W$ is the height of the free energy wall, $\epsilon^2$ is the gradient energy coefficient, and $\gamma$ is a phenomenological parameter that quantifies the free energy increase at overlapping interfaces.
\section{Governing Equations}
The evolution of the system is determined via the constrained minimization of the free energy with respect to $\eta$  and  $\psi$ coupled to the electrochemical kinetics at the interface and the transport of ionic species in the electrolyte. The order parameter is evolved according to the advective Cahn-Hilliard equation:
\begin{equation}
\label{CH_eta}
 \frac{\partial \eta}{\partial t} = \nabla \cdot \left[ M(\psi) \nabla \frac{\delta \Pi}{\delta \eta}\right] + v |\nabla \psi|,
 \end{equation}
\noindent where $M$ is the Cahn-Hilliard mobility coefficient and $v$ is the velocity of the interface normal to the surface, which is related to the reaction current at the interface via Faraday's law of electrolysis:
\begin{equation}
\label{int_vel}
v=-\frac{V_M i_{rxn}}{z_M F}
 \end{equation}
\noindent where $V_M$ is the molar volume of the metal, $i_{rxn}$ is the reaction current density, $z_M$ is the dissolved metal cation charge number, and $F$ is Faraday's constant. The interfacial mobility, $M$, is set as a function of the local current density,
\begin{equation}
M=2\frac{V_M |i_{rxn}| \psi}{z_M F}\sqrt{\frac{2 \epsilon^2}{W}}
\end{equation}
The evolution of $\psi$ employs the advective Cahn-Hilliard equation, but with the opposite sign of the advective source term in Eq. \eqref{CH_eta}:
\begin{equation}
\label{CH_psi}
 \frac{\partial \psi}{\partial t} = \nabla \cdot \left[ M(\psi) \nabla \frac{\delta \Pi}{\delta \psi} \right] - v |\nabla \psi|.
 \end{equation}
 The smoothed boundary method (SBM)  [H.-C. Yu, H.-Y. Chen, and K. Thornton, Model. Simul. Mater. Sci. Eng. {\bf 20}, 075008
(2012)] is employed to confine the transport of ionic species to the electrolyte region, where $\psi > 0$. Thus, the ionic transport equations are coupled to the phase field method by using  $\psi$ as the domain parameter. The electrolyte is assumed to be composed of three ionic species: an effective cation species for the metal (M) dissolved at the interface, a supporting cation (+), and a supporting anion (-). The electrolyte is assumed to be electroneutral, which implies that that the following constraint must be satisfied,
\begin{equation}
\label{electroneutrality}
z_M c_M + z_+ c_+ + z_-c_- = 0.
 \end{equation}
Thus, the concentrations are not independent of each other only two additional equations are required to describe the concentration evolution for all three species.  The concentrations of the effective metal cation and the supporting electrolyte cation are directly solved and the supporting anion concentration is eliminated via Eq. \eqref{electroneutrality}. The SBM reformulated governing equation for each species, $i$, which includes diffusion and migration effects, is given by
\begin{equation}
\label{conc_dynamics}
\frac{\partial c_i}{\partial t}=\frac{1}{\psi} \nabla \cdot (\psi D_i \nabla c_i) + \frac{1}{\psi} \left( \frac{z_i F}{RT} \nabla \cdot (\psi D_i c_i \nabla \Phi) \right)
+ \frac{|\nabla \psi|}{\psi} \left( \frac{i_{rxn}}{z_i F} \right)
\end{equation}
\noindent where $c_i =c_M, c_+$ and $\Phi$ is the electrostatic potential.  The last term of the right hand side of Eq. \eqref{conc_dynamics} is zero for $c_i =c_+$  because $i_{rxn}$ is the reaction current density for the dissolution of the metal (M), which does not involve the supporting cation (+) and anion (-). The SBM reformulated governing equation for the potential is
\begin{equation}
\nabla \cdot (\psi \kappa \nabla \Phi) = F \nabla \cdot \left[ \psi \left( z_M (D_- - D_M) \nabla c_M + z_+ (D_- - D_+)  \nabla c_+ \right) \right] - |\nabla \psi | i_{rxn},
\end{equation}
where 
\begin{equation}
\kappa =\frac{F^2}{RT}\left[ z_Mc_M\left( z_MD_M-z_-D_-\right) + z_+c_+\left( z_+D_+-z_-D_-\right) \right].
\end{equation}
The reaction current density is obtained via a Butler-Volmer type kinetic expression that includes a maximum current density that accounts for both how far the electrolyte is from saturation at a point in time as well as the rate of transport of ions into the electrolyte:
\begin{equation}
\frac{i_{rxn}}{i_{corr}}= \left( 1 - \frac{i_{rxn}}{i_{max,c}} \right) \exp \left( \frac{z_M (1-\beta) F}{RT}\xi \right),
\end{equation}
where $i_{max,c}$ is the maximum reaction current density, $\beta$ is the the charge transfer symmetry factor and $\xi$ is the overpotential, defined as $\xi=V_s-E_{corr}-\Phi$. The parameters $i_{corr}$ and $E_{corr}$ are the corrosion current density and corrosion potential, respectively, and $V_s$ is the applied potential. The maximum reaction current density is given by
\begin{equation}
i_{max,c}=\left( \frac{z_M F}{1-c_M V_M} \right) \left[ \frac{2\delta}{\tau}(c_{M,sat} - c_M ) + \left( D_M\nabla c_M +z_M\frac{F}{RT}D_M c_M \nabla \Phi \right) \cdot \mathbf{n} \right], 
\end{equation}
where $c_{M,sat}$ is the saturation concentration of the metal ions in solution and $2\delta/\tau$ is a characteristic velocity of ion transport across the diffuse interface. This velocity is given by the largest value between the characteristic velocities of diffusion and of migration:
\begin{equation}
\frac{2\delta}{\tau}=\max \left( \frac{D_M}{2\delta}, \left| \frac{z_M D_M F \nabla \Phi \cdot \mathbf{n}}{RT} \right| \right).
\end{equation}
\section{Time Discretization}
The fields $\eta$, $\psi$, $c_M$ and $c_+$ are solved using an explicit Euler method for time integration. For  $\eta$ and $\psi$ we employ a splitting strategy to transform Eqs. \eqref{CH_eta} and \eqref{CH_psi} into second order PDEs:
\begin{equation}
\label{eta_eq}
 \frac{\partial \eta}{\partial t} = \nabla \cdot \left( M \nabla  \mu_\eta \right) + v |\nabla \psi|
\end{equation}
and
\begin{equation}
\label{psi_eq}
\frac{\partial \psi}{\partial t} = \nabla \cdot \left( M \nabla \mu_\psi \right) - v |\nabla \psi|,
\end{equation}
where 
\begin{equation}
\label{mueta_eq}
\mu_\eta=W \left( \eta^3 - \eta +2 \gamma \eta \psi^2 \right) -\epsilon^2 \nabla^2 \eta
\end{equation}
and
\begin{equation}
\label{mupsi_eq}
\mu_\psi=W \left( \psi^3 - \psi +2 \gamma \psi  \eta^2 \right) -\epsilon^2 \nabla^2 \psi.
\end{equation}
Considering forward Euler explicit time stepping, the time-discretized version of the equations above are
\begin{equation}
\label{eta_eq_td}
\eta^{n+1} = \eta^n + \Delta t \left[ \nabla \cdot \left( M \nabla \mu_\eta^n \right) + v |\nabla \psi^n| \right],
\end{equation}
\begin{equation}
\label{psi_eq_td}
\psi^{n+1} = \psi^n +  \Delta t \left[  \nabla \cdot \left( M \nabla \mu_\psi^n \right) - v |\nabla \psi^n | \right],
\end{equation}
\begin{equation}
\label{mueta_eq_td}
\mu_\eta^{n+1}=f_\eta^n -\epsilon^2 \nabla^2 \eta^n,
\end{equation}
and
\begin{equation}
\label{mupsi_eq_td}
\mu_\psi^{n+1}=f_\psi^n -\epsilon^2 \nabla^2 \psi^n,
\end{equation}
where $f_\eta=W \left( \eta^3 - \eta +2 \gamma \eta \psi^2 \right)$ and $f_\psi= W \left( \psi^3 - \psi +2 \gamma \psi  \eta^2 \right)$.\\
The discretized equation for the ion concentrations are
\begin{equation}
\label{conc_eq_td}
c_i ^{n+1}= c_i^n  + \Delta t  \left[ \frac{1}{\psi^n} \nabla \cdot (\psi^n D_i \nabla c_i^n) + \frac{1}{\psi^n} \left( \frac{z_i F}{RT} \nabla \cdot (\psi^n D_i c_i^n \nabla \Phi^n) \right)+ \frac{|\nabla \psi^n|}{\psi^n} \left( \frac{i_{rxn}}{z_i F} \right) \right].
\end{equation}
The electrostatic potential is assumed to be in equilibrium throughout the simulation and needs to be solved as a non-linear time-independent equation.
\section{Weak Formulation}
For the weak formulation of time-discretized equations \eqref{eta_eq_td}-\eqref{conc_eq_td} only RHS terms need to be specified (see \href{https://prisms-center.github.io/phaseField/doxygen_files/app_files.html}{\textcolor{blue}{User Manual}} for details)
\begin{equation}
\int_{\Omega} \omega \eta^{n+1} dV = \int_{\Omega} \omega \underbrace{\left( \eta^n +\Delta t v  |\psi^n|\right)}_{r_\eta}dV + \int_{\Omega} \nabla \omega \cdot  \underbrace{\left(-  \Delta t M \nabla \mu_\eta^n \right)}_{r_{\eta x}}dV
\end{equation}'
\begin{equation}
\int_{\Omega} \omega \mu_\eta^{n+1} dV = \int_{\Omega} \omega \underbrace{f_\eta^n}_{r_{\mu \eta}}dV + \int_{\Omega} \nabla \omega \cdot  \underbrace{(\epsilon^2\nabla\eta^n)}_{r_{\mu \eta x}}dV
\end{equation}
\begin{equation}
\int_{\Omega} \omega \psi^{n+1} dV = \int_{\Omega} \omega \underbrace{\left( \psi^n -\Delta t v  |\psi^n|\right)}_{r_\psi} dV+ \int_{\Omega} \nabla \omega \cdot  \underbrace{\left(-  \Delta t M \nabla \mu_\psi^n \right)}_{r_{\psi x}}dV
\end{equation}
\begin{equation}
\int_{\Omega} \omega \mu_\psi^{n+1} dV = \int_{\Omega} \omega \underbrace{f_\psi^n}_{r_{\mu \psi}}dV + \int_{\Omega} \nabla \omega \cdot  \underbrace{(\epsilon^2\nabla\psi^n)}_{r_{\mu \psi x}}dV
\end{equation}
\begin{align}
\int_{\Omega} \omega c_i^{n+1} dV &=\int_{\Omega} \omega \underbrace{\left( c_i^n+\frac{\Delta t D_i}{\psi^n}\nabla\psi^n\cdot\nabla c_i^n 
+\frac{\Delta t D_i z_i F}{RT\psi^n}\nabla\psi^n\cdot(c_i^n\nabla \Phi^n) +\frac{\Delta t}{z_i F\psi^n} |\nabla\psi^n| i_{rxn}\right)}_{r_{ci}} dV\\ 
&+ \int_{\Omega} \nabla \omega \cdot  \underbrace{\left( -\Delta t D_i \nabla c_i^n -\frac{\Delta t D_i z_i F}{RT} c_i^n\nabla \Phi^n \right)}_{r_{cix}}dV
\end{align}

For the time-independent electrostatic potential, we need to specify LHS and RHS terms:
\begin{align}
&\int_{\Omega} \omega  \underbrace{\left( \frac{\partial i_{rxn}}{\partial \Phi}^n |\nabla\psi^n| \Delta\Phi\right)}_{r_{\Delta \Phi}} dV 
+ \int_{\Omega} \nabla\omega \cdot  \underbrace{\left( -\psi^n \kappa  \nabla(\Delta\Phi)\right)}_{r_{\Delta \Phi x}} dV =\\
&\int_{\Omega} \omega  \underbrace{\left( -|\nabla\psi^n| i_{rxn} \right)}_{r_\Phi} dV + \int_{\Omega} \nabla\omega \cdot  \underbrace{\left( -F  \left[ \left( z_M (D_- - D_M) \nabla c_M^n + z_+ (D_- - D_+)  \nabla c_+^n \right) \right] +\psi^n \kappa \nabla\Phi^n \right)}_{r_\Phi x} dV
\end{align}

The above values of are used to define the equation terms in the input file $applications/corrosion/equations.cc$.
