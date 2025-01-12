# KKS Phase Field Model of Precipitate Evolution coupled with nucleation (October 2, 2024)

The Nucleation Model application for PRISMS-PF  incorporates a stochastic method to add nuclei to the KKS phase field model for precipitate growth. Nuclei are seeded throughout the time evolution of the system based on a probability that depends on the local solute supersaturation. This document is divided in two sections. In the first section, the phase field model formulation for precipitate evolution in a binary alloy (without elastic effects) is presented. In the second section the nucleation method is presented. 

\section{Precipitate Evolution}
\subsection{Variational formulation}
In the absence of elastic effects total free energy of the 2-component system (neglecting boundary terms) is of the form,
\begin{equation}
\Pi(c, \eta) = \int_{\Omega} f(c, \eta) ~dV 
\end{equation}
where $c$ is the concentration of the $\beta$ phase and  $\eta$ is the set of structural order parameters. The free energy density, $f$, is given by
\begin{equation}
 f(c, \eta) =   f_{chem}(c, \eta) + f_{grad}(\eta)
\end{equation}
where
\begin{equation}
f_{chem}(c, \eta) = f_{\alpha}(c,\eta) \left( 1- H(\eta)\right) + f_{\beta}(c,\eta) H(\eta)+ W f_{Landau}(\eta)
\end{equation}
and
\begin{equation}
f_{grad}(\eta) = \frac{1}{2} \kappa | \nabla \eta |^2 \\
\end{equation}

In the KKS model (Kim 1999), the interfacial region is modeled as a mixture of the $\alpha$ and $\beta$ phases with concentrations $c_{\alpha}$ and $c_{\beta}$, respectively. The homogenous free energies for each phase, $f_{\alpha}$ and $f_{\beta}$ in this case, are typically given as functions of $c_{\alpha}$ and $c_{\beta}$, rather than directly as functions of $c$ and $\eta_p$. Thus, $f_{chem}(c, \eta)$ can be rewritten as 
\begin{equation}
f_{chem}(c, \eta) = f_{\alpha}(c_\alpha) \left( 1- H(\eta)\right) + f_{\beta}(c_\beta) H(\eta)+ W f_{Landau}(\eta)
\end{equation}

The concentration in each phase is determined by the following system of equations:
\begin{gather}
c =  c_{\alpha} \left( 1- H(\eta)\right) + c_{\beta} H(\eta) \\
\frac{\partial f_{\alpha}(c_{\alpha})}{\partial c_{\alpha}} = \frac{\partial f_{\beta}(c_{\beta})}{\partial c_{\beta}}
\end{gather}

Given the following parabolic functions for the single-phase homogenous free energies:
\begin{gather}
f_{\alpha}(c_{\alpha}) = A_{2} c_{\alpha}^2 + A_{1} c_{\alpha} + A_{0} \\
f_{\beta}(c_{\beta}) = B_{2} c_{\beta}^2 + B_{1} c_{\beta} + B_{0}
\end{gather}
the single-phase concentrations are:
\begin{gather}
c_{\alpha} = \frac{ B_2 c + \frac{1}{2} (B_1 - A_1) H(\eta) }{A_2 H(\eta) + B_2 \left( 1- H(\eta)\right) } \\
c_{\beta} =  \frac{ A_2 c + \frac{1}{2} (A_1 - B_1) \left[1-H(\eta)\right] }{A_2  H(\eta) + B_2 \left[ 1- H(\eta)\right] } 
\end{gather}

\subsection{Required inputs}
\begin{itemize}
\item $f_{\alpha}(c_{\alpha}), f_{\beta}(c_{\beta})$ - Homogeneous chemical free energy of the components of the binary system, example form given above
\item $f_{Landau}(\eta)$ - Landau free energy term that controls the interfacial energy. Example form given in Appendix I
\item $W$ - Barrier height for the Landau free energy term, used to control the thickness of the interface 
\item $H(\eta)$ - Interpolation function for connecting the $\alpha$ phase and the $\beta$ phase. Example form given in Appendix I
\item $\Bkappa^{\eta_p}$  - gradient penalty coefficient for the $\alpha - \beta$ interface
\end{itemize}
In addition, to drive the kinetics, we need:
\begin{itemize}
\item $M$  - mobility value for the concentration field
\item $L$  - mobility value for the structural order parameter field
\end{itemize}

\subsection{Variational treatment}
We obtain chemical potentials for the concentration and the structural order parameter by taking variational derivatives of $\Pi$:
\begin{align}
  \mu_{c}  &= f_{\alpha,c} \left( 1- H(\eta)\right) +f_{\beta,c} H(\eta) \\
  \mu_{\eta}  &= \left[ f_{\beta}-f_{\alpha} -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}} \right] H(\eta)_{,\eta} + W f_{Landau,\eta}- \kappa\nabla^2\eta
\end{align}

\subsection{Kinetics}
Now the PDE for Cahn-Hilliard dynamics is given by:
\begin{align}
  \frac{\partial c}{\partial t} &= ~\grad \cdot \left( \frac{1}{f_{,cc}}M \grad \mu_c \right) \label{CH_eqn}
  \end{align}
  where $M$ is a constant mobility and the factor of $\frac{1}{f_{,cc}}$ is added to guarentee constant diffusivity in the two phases. The PDE for Allen-Cahn dynamics is given by:
  \begin{align}
    \frac{\partial \eta}{\partial t} &= - L \mu_{\eta_p} \label{AC_eqn}
\end{align}
where  $L$ is a constant mobility. 

\subsection{Time discretization}
Using forward Euler explicit time stepping, equations \ref{CH_eqn} and \ref{AC_eqn} become:
\begin{align}
c^{n+1} = c^{n}+\Delta t \left[\nabla \cdot \left(\frac{1}{f_{,cc}} M \nabla \mu_c \right) \right]\\
\eta_p^{n+1} = \eta_p^n -\Delta t L \mu_{\eta_p}
\end{align}

\subsection{Weak formulation}
Writing equations \ref{CH_eqn} and \ref{AC_eqn} in the weak form, with the arbirary variation given by $w$ yields:
\begin{align}
\int_\Omega w c^{n+1} dV &= \int_\Omega wc^{n}+w  \Delta t \left[\nabla \cdot \left(\frac{1}{f_{,cc}}  M \nabla \mu_c \right) \right] dV \label{CH_weak} \\
%&= \int_\Omega w\underbrace{c^{n}}_{r_c}+\nabla w \cdot (\underbrace{\Delta t  M \nabla \mu_c}_{r_{cx}} ) dV \\
\int_\Omega w \eta^{n+1} dV &= \int_\Omega w \eta^{n}-w  \Delta t L \mu_{\eta} dV  \label{AC_weak}
%&= \int_\Omega w\underbrace{c^{n}}_{r_c}+\nabla w \cdot (\underbrace{\Delta t  M \nabla \mu_c}_{r_{cx}} ) dV 
\end{align}

The expression of $\frac{1}{f_{,cc}} \mu_c$ can be written as:
\begin{equation}
\frac{1}{f_{,cc}}  \nabla \mu_c =  \nabla c + (c_{\alpha}-c_{\beta}) H(\eta)_{,\eta} \nabla \eta  \\
\end{equation}

Applying the divergence theorem to equation \ref{CH_weak}, one can derive the residual terms $r_c$ and $r_{cx}$:
\begin{equation}
\int_\Omega w c^{n+1} dV = \int_\Omega w\underbrace{c^{n}}_{r_c}+\nabla w \cdot (\underbrace{-\Delta t  M \frac{1}{f_{,cc}} \nabla \mu_c}_{r_{cx}} ) dV
\end{equation}

Expanding $\mu_{\eta}$ in equation \ref{AC_weak} and applying the divergence theorem yields the residual terms $r_{\eta}$ and $r_{\eta x}$:
\begin{equation}
\begin{split}
\int_\Omega w \eta^{n+1} dV = &\int_\Omega w \Bigg\{\underbrace{\eta^{n}-\Delta t L \bigg[(f_{\beta}-f_{\alpha})H(\eta^n)_{,\eta} -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}}H(\eta^n)_{,\eta} + W f_{Landau,\eta}}_{r_{\eta}}\\ 
&+ \nabla w \cdot (\underbrace{-\Delta t  L \kappa \nabla \eta^n}_{r_{\eta x}} ) dV 
\end{split}
\end{equation}

\section{Nucleation method}

We follow the same approach as Jokisaari and Thornton [Comput. Mater. Sci. {\bf 112}, 128-138 (2016)] which consists of adding nuclei throughout a phase field simulation based on a probability that depends on the local supersaturation. This probability is calculated every fixed number of time steps and for every element of the computational domain. In each nucleation event, nucleation is triggered at a point within the $\alpha$ phase. Each nucleus is then added to the system by modifying the order parameter to it's $\beta$ phase value within a small domain around the selected nucleation center. This domain can be spherical/circular or ellipsoidal/elliptical.  

\subsection{Nucleation rate}

From classical nucleation theory, the nucleation rate for critical nuclei $J^*$ is given by
\\
\begin{equation}
\label{nuc_rate_full}
J^*(\mathbf{r},t)=Zn\beta^*\exp \left( -\frac{\Delta G^*}{k_B T} \right) \exp \left( -\frac{\tau}{t} \right),
\end{equation}
\\
where $Z$ is the Zeldovich factor, $n$ is the number of nucleation sites per volume, $\beta^*$ is the frequency at which a critical nucleus becomes supercritical, $\Delta G^*$ is the nucleation energy barrier, $k_B$ is the Boltzmann constant, $T$ is the temperature, $t$ is time and $\tau$ is the incubation time. It can be shown that, in the dilute limit and for constant temperature, Eq.~\eqref{nuc_rate_full} can be simplified by grouping approximately constant terms in both the exponential and pre-exponential factors:
\\
\begin{equation}
\label{nuc_rate_simp}
J^*(\mathbf{r},t)=k_1\exp \left( -\frac{k_2}{(\Delta c)^{d-1}} \right) \exp \left(-\frac{\tau}{t} \right),
\end{equation}
\\
where  $k_1$ and $k_2$ are now taken as constant parameters, $\Delta c=c(\mathbf{r},t)-c_\alpha^{eq}$ is the local supersaturation in the $\alpha$ phase and $d$ is the dimensionality of the system ({\it e.g.} $d=2$ or $d=3$).\\

\subsection{Nucleation probability}

Considering  $J^*$ to be approximately constant within a small volume, $\Delta V$, and for a small time interval, $\Delta t$, the probability that at least one nucleation event occurs in $\Delta V$ within $\Delta t$ is given by [Simmons et al., Scripta Mater. {\bf 43}, 935 (2000)]
\\
\begin{equation}
\label{nuc_prob}
P(\mathbf{r},t) = 1 - \exp \left( -J^* \Delta V \Delta t \right)
\end{equation}
\subsection{Hold time}

After each nucleus is added, there is a `hold' time interval, $\Delta t_h$, during which the order parameter value is fixed within a small window that encompasses the new nucleus. The purpose of this hold time is to allow the concentration to evolve within the nucleus to a value close to the coexistance composition for $\beta$ phase, and therefore, to create small a solute depleted zone around the nucleus. After the hold time, the nucleus is allowed to evolve into a precipitate.

\subsection{Required nucleation inputs}
\begin{itemize}
\item $k_1$ - Constant pre-exponential factor in Eq.~\eqref{nuc_rate_simp}
\item $k_2$ - Parameter that groups all constant terms of the first exponential factor in Eq.~\eqref{nuc_rate_simp} 
\item $\tau$ - Incubation time constant in Eq.~\eqref{nuc_rate_simp}
\item $\Delta t_h$ - Nucleation hold time.
\end{itemize}
Dimensions (ellipsoidal semiaxes) of precipitate seeds
\begin{itemize}
\item a - semiaxis in the x-direction
\item b - semiaxis in the y-direction
\item c - semiaxis in the x-direction
\end{itemize}

\section*{Appendix I: Example functions for $f_{\alpha}$, $f_{\beta}$, $f_{Landau}$, $H(\eta)$ }
\begin{gather}
f_{\alpha}(c_{\alpha}) = A_{2} c_{\alpha}^2 + A_{1} c_{\alpha} + A_{0} \\
f_{\beta}(c_{\beta}) = B_{2} c_{\beta}^2 + B_{1} c_{\beta} + B_{0} \\
f_{Landau}(\eta) = \eta^2  - 2\eta^3 +  \eta^4\\
H(\eta) = 3 \eta^2 - 2 \eta^3
\end{gather}
