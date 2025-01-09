## KKS Phase Field Model of Precipitate Evolution (September 14th, 2023)
 
The model employed in this application is described in detail in the article:
DeWitt et al., Misfit-driven $\beta'''$ precipitate composition and morphology in Mg-Nd alloys,
Acta Materialia **137**, 378-389 (2017).
 
## Variational formulation
The total free energy of the system (neglecting boundary terms) is of the form,

$$
\begin{equation}
\Pi(c, \eta_1, \eta_2, \eta_3, \epsilon) = \int_{\Omega} f(c, \eta_1, \eta_2, \eta_3, \epsilon) ~dV 
\end{equation}
$$

where $c$ is the concentration of the $\beta$ phase, $\eta_p$ are the structural order parameters and $\varepsilon$ is the small strain tensor. $f$, the free energy density is given by

$$
\begin{equation}
 f(c, \eta_1, \eta_2, \eta_3, \epsilon) =   f_{chem}(c, \eta_1, \eta_2, \eta_3) + f_{grad}(\eta_1, \eta_2, \eta_3) + f_{elastic}(c,\eta_1, \eta_2, \eta_3,\epsilon)
\end{equation}
$$

where

$$
\begin{equation}
f_{chem}(c, \eta_1, \eta_2, \eta_3) = f_{\alpha}(c,\eta_1, \eta_2, \eta_3) \left( 1- \sum_{p=1}^3 H(\eta_p)\right) + f_{\beta}(c,\eta_1, \eta_2, \eta_3) \sum_{p=1}^3 H(\eta_p)+ W f_{Landau}(\eta_1, \eta_2, \eta_3)
\end{equation}
$$

$$
\begin{equation}
f_{grad}(\eta_1, \eta_2, \eta_3) = \frac{1}{2} \sum_{p=1}^3 \Bkappa^{\eta_p}_{ij} \eta_{p,i}  \eta_{p,j}
\end{equation}
$$

$$
\begin{gather}
f_{elastic}(c,\eta_1, \eta_2, \eta_3,\epsilon) = \frac{1}{2} \bC_{ijkl}(\eta_1, \eta_2, \eta_3)  \left( \varepsilon_{ij} - \varepsilon ^0_{ij}(c, \eta_1, \eta_2, \eta_3) \right)\left( \varepsilon_{kl} - \varepsilon^0_{kl}(c, \eta_1, \eta_2, \eta_3)\right) 
\end{gather}
$$

$$
\begin{gather}
\varepsilon^0(c, \eta_1, \eta_2, \eta_3) = H(\eta_1) \varepsilon^0_{\eta_1} (c_{\beta})+ H(\eta_2) \varepsilon^0_{\eta_2} (c_{\beta}) + H(\eta_3) \varepsilon^0_{\eta_3} (c_{\beta}) 
\bC(\eta_1, \eta_2, \eta_3) = H(\eta_1) \bC_{\eta_1}+ H(\eta_2) \bC_{\eta_2} + H(\eta_3) \bC_{\eta_3} + \left( 1- H(\eta_1)-H(\eta_2)-H(\eta_3)\right)  \bC_{\alpha}
\end{gather}
$$

Here $\varepsilon^0_{\eta_p}$ are the composition dependent stress free strain transformation tensor corresponding to each structural order parameter, which is a function of the $\beta$ phase concentration, $c_{\beta}$, defined below.

In the KKS model (Kim 1999), the interfacial region is modeled as a mixture of the $\alpha$ and $\beta$ phases with concentrations $c_{alpha}$ and $c_{beta}$, respectively. The homogenous free energies for each phase, $f_{\alpha}$ and $f_{\beta}$ in this case, are typically given as functions of $c_{\alpha}$ and $c_{\beta}$, rather than directly as functions of $c$ and $\eta_p$. Thus, $f_{chem}(c, \eta_1, \eta_2, \eta_3)$ can be rewritten as 
\begin{equation}
f_{chem}(c, \eta_1, \eta_2, \eta_3) = f_{\alpha}(c_{\alpha}) \left( 1- \sum_{p=1}^3 H(\eta_p)\right) + f_{\beta}(c_{\beta}) \sum_{p=1}^3 H(\eta_p)+ W f_{Landau}(\eta_1, \eta_2, \eta_3) 
\end{equation}

The concentration in each phase is determined by the following system of equations:
\begin{gather}
c =  c_{\alpha} \left( 1- \sum_{p=1}^3 H(\eta_p)\right) + c_{\beta} \sum_{p=1}^3 H(\eta_p) \\
\frac{\partial f_{\alpha}(c_{\alpha})}{\partial c_{\alpha}} = \frac{\partial f_{\beta}(c_{\beta})}{\partial c_{\beta}}
\end{gather}

Given the following parabolic functions for the single-phase homogenous free energies:
\begin{gather}
f_{\alpha}(c_{\alpha}) = A_{2} c_{\alpha}^2 + A_{1} c_{\alpha} + A_{0} \\
f_{\beta}(c_{\beta}) = B_{2} c_{\beta}^2 + B_{1} c_{\beta} + B_{0}
\end{gather}
the single-phase concentrations are:
\begin{gather}
c_{\alpha} = \frac{ B_2 c + \frac{1}{2} (B_1 - A_1) \sum_{p=1}^3 H(\eta_p) }{A_2 \sum_{p=1}^3 H(\eta_p) + B_2 \left( 1- \sum_{p=1}^3 H(\eta_p)\right) } \\
c_{\beta} =  \frac{ A_2 c + \frac{1}{2} (A_1 - B_1) \left[1-\sum_{p=1}^3 H(\eta_p)\right] }{A_2 \sum_{p=1}^3 H(\eta_p) + B_2 \left[ 1- \sum_{p=1}^3 H(\eta_p)\right] } 
\end{gather}

\section{Required inputs}
\begin{itemize}
\item $f_{\alpha}(c_{\alpha}), f_{\beta}(c_{\beta})$ - Homogeneous chemical free energy of the components of the binary system, example form given above
\item $f_{Landau}(\eta_1, \eta_2, \eta_3)$ - Landau free energy term that controls the interfacial energy and prevents precipitates with different orientation varients from overlapping, example form given in Appendix I
\item $W$ - Barrier height for the Landau free energy term, used to control the thickness of the interface 
\item $H(\eta_p)$ - Interpolation function for connecting the $\alpha$ phase and the $p^{th}$ orientation variant of the $\beta$ phase, example form given in Appendix I
\item $\Bkappa^{\eta_p}$  - gradient penalty tensor for the $p^{th}$ orientation variant of the $\beta$ phase
\item $\bC_{\eta_p}$ - fourth order elasticity tensor (or its equivalent second order Voigt representation) for the $p^{th}$ orientation variant of the $\beta$ phase
\item $\bC_{\alpha}$ - fourth order elasticity tensor (or its equivalent second order Voigt representation) for the $\alpha$ phase
\item $\varepsilon^0_{\eta_p}$ - stress free strain transformation tensor for the $p^{th}$ orientation variant of the $\beta$ phase
\end{itemize}
In addition, to drive the kinetics, we need:
\begin{itemize}
\item $M$  - mobility value for the concentration field
\item $L$  - mobility value for the structural order parameter field
\end{itemize}

\section{Variational treatment}
%From the variational derivatives given in Appendix II, we obtain the chemical potentials for the concentration and the structural order parameters:
We obtain chemical potentials for the chemical potentials for the concentration and the structural order parameters by taking variational derivatives of $\Pi$:
\begin{align}
  \mu_{c}  &= f_{\alpha,c} \left( 1- H(\eta_1)-H(\eta_2)-H(\eta_3)\right) +f_{\beta,c} \left(  H(\eta_1)  + H(\eta_2) + H(\eta_3) \right)  + \bC_{ijkl} (- \varepsilon^0_{ij,c}) \left( \varepsilon_{kl} - \varepsilon^0_{kl}\right) \\
  \mu_{\eta_p}  &= [ f_{\beta}-f_{\alpha} -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}}] H(\eta_p)_{,\eta_p} + W f_{Landau,\eta_p}- \Bkappa^{\eta_p}_{ij} \eta_{p,ij} + \bC_{ijkl} (- \varepsilon^0_{ij,\eta_p}) \left( \varepsilon_{kl} - \varepsilon^0_{kl}\right) + \frac{1}{2} \bC_{ijkl,\eta_p} \left( \varepsilon_{ij} - \varepsilon ^0_{ij} \right) \left( \varepsilon_{kl} - \varepsilon^0_{kl}\right)
\end{align}

\section{Kinetics}
Now the PDE for Cahn-Hilliard dynamics is given by:
\begin{align}
  \frac{\partial c}{\partial t} &= ~\grad \cdot \left( \frac{1}{f_{,cc}}M \grad \mu_c \right) \label{CH_eqn}
  \end{align}
  where $M$ is a constant mobility and the factor of $\frac{1}{f_{,cc}}$ is added to guarentee constant diffusivity in the two phases. The PDE for Allen-Cahn dynamics is given by:
  \begin{align}
    \frac{\partial \eta_p}{\partial t} &= - L \mu_{\eta_p} \label{AC_eqn}
\end{align}
where  $L$ is a constant mobility. 

\section{Mechanics}
Considering variations on the displacement $u$ of the from $u+\epsilon w$, we have
\begin{align}
\delta_u \Pi &=  \int_{\Omega}   \grad w :  \bC(\eta_1, \eta_2, \eta_3) : \left( \varepsilon - \varepsilon^0(c,\eta_1, \eta_2, \eta_3)\right) ~dV = 0 \\
\end{align}
where $\Bsigma = \bC(\eta_1, \eta_2, \eta_3) : \left( \varepsilon - \varepsilon^0(c,\eta_1, \eta_2, \eta_3)\right)$ is the stress tensor. \\

\section{Time discretization}
Using forward Euler explicit time stepping, equations \ref{CH_eqn} and \ref{AC_eqn} become:
\begin{align}
c^{n+1} = c^{n}+\Delta t \left[\nabla \cdot \left(\frac{1}{f_{,cc}} M \nabla \mu_c \right) \right]\\
\eta_p^{n+1} = \eta_p^n -\Delta t L \mu_{\eta_p}
\end{align}

\section{Weak formulation}
Writing equations \ref{CH_eqn} and \ref{AC_eqn} in the weak form, with the arbirary variation given by $w$ yields:
\begin{align}
\int_\Omega w c^{n+1} dV &= \int_\Omega wc^{n}+w  \Delta t \left[\nabla \cdot \left(\frac{1}{f_{,cc}}  M \nabla \mu_c \right) \right] dV \label{CH_weak} \\
%&= \int_\Omega w\underbrace{c^{n}}_{r_c}+\nabla w \cdot (\underbrace{\Delta t  M \nabla \mu_c}_{r_{cx}} ) dV \\
\int_\Omega w \eta_p^{n+1} dV &= \int_\Omega w \eta_p^{n}-w  \Delta t L \mu_{\eta_p} dV  \label{AC_weak}
%&= \int_\Omega w\underbrace{c^{n}}_{r_c}+\nabla w \cdot (\underbrace{\Delta t  M \nabla \mu_c}_{r_{cx}} ) dV 
\end{align}

The expression of $\frac{1}{f_{,cc}} \mu_c$ can be written as:
\begin{equation}
\begin{split}
\frac{1}{f_{,cc}}  \nabla \mu_c = & \nabla c + (c_{\alpha}-c_{\beta}) \sum_{p=1}^3 H(\eta_p)_{,\eta_p} \nabla \eta_p  \\
&+ \frac{1}{f_{,cc}} \left[ \sum_{p=1}^3 (C_{ijkl}^{\eta_p} - C_{ijkl}^{\alpha} )\nabla \eta_p H(\eta_p)_{,\eta_p} \right](-\epsilon_{ij,c}^0)(\epsilon_{ij} - \epsilon_{ij}^0) \\
&- \frac{1}{f_{,cc}} C_{ijkl} \left[  \sum_{p=1}^3 \left( H(\eta_p)_{,\eta_p} \epsilon_{ij,c}^{0\eta_p} + \sum_{q=1}^3 \left( H(\eta_p) \epsilon_{ij,c\eta_q}^{0\eta_p} \right) \right) \nabla \eta_p + H(\eta_p) \epsilon_{ij,cc}^{0\eta_p} \nabla c \right](\epsilon_{kl}-\epsilon_{kl}^0)\\
&+ \frac{1}{f_{,cc}} C_{ijkl} (-\epsilon_{ij,c}^0) \left[ \nabla \epsilon_{kl} -  \left( \sum_{p=1}^3 \left(H(\eta_p)_{,\eta_p} \epsilon_{kl}^{0\eta_p} -\sum_{q=1}^3 \epsilon_{kl,\eta_q}^{\eta_q} H(\eta_q) \right)\nabla \eta_p + H(\eta_p) \epsilon_{kl,c}^{0\eta_p} \nabla c \right) \right]
\end{split}
\end{equation}

Applying the divergence theorem to equation \ref{CH_weak}, one can derive the residual terms $r_c$ and $r_{cx}$:
\begin{equation}
\int_\Omega w c^{n+1} dV = \int_\Omega w\underbrace{c^{n}}_{r_c}+\nabla w \cdot (\underbrace{-\Delta t  M \frac{1}{f_{,cc}} \nabla \mu_c}_{r_{cx}} ) dV
\end{equation}

Expanding $\mu_{\eta_p}$ in equation \ref{AC_weak} and applying the divergence theorem yields the residual terms $r_{\eta_p}$ and $r_{\eta_p x}$:
\begin{equation}
\begin{split}
\int_\Omega w \eta_p^{n+1} dV &= \\
&\int_\Omega w \Bigg\{\underbrace{\eta_p^{n}-\Delta t L \bigg[(f_{\beta}-f_{\alpha})H(\eta_p^n)_{,\eta_p} -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}}H(\eta_p^n)_{,\eta_p} + W f_{Landau,\eta_p}}_{r_{\eta_p}}  \\ 
&\underbrace{ -C_{ijkl} \left( H(\eta_p)_{,\eta_p} \epsilon_{ij}^{0 \eta_p}\right)\left(\epsilon_{kl} - \epsilon_{kl}^{0} \right) + \frac{1}{2} \left[ (C_{ijkl}^{\eta_p} - C_{ijkl}^{\alpha}) H(\eta_p)_{,\eta_p} \right] \left(\epsilon_{ij} - \epsilon_{ij}^{0} \right) \left(\epsilon_{kl} - \epsilon_{kl}^{0} \right) \bigg] }_{r_{\eta_p}~cont.} \Bigg\} \\
&+ \nabla w \cdot (\underbrace{-\Delta t  L \Bkappa^{\eta_p}_{ij} \eta_{p,i}^n}_{r_{\eta_p x}} ) dV 
\end{split}
\end{equation}

\section{Appendix I: Example functions for $f_{\alpha}$, $f_{\beta}$, $f_{Landau}$, $H(\eta_p)$ }
\begin{gather}
f_{\alpha}(c_{\alpha}) = A_{2} c_{\alpha}^2 + A_{1} c_{\alpha} + A_{0} \\
f_{\beta}(c_{\beta}) = B_{2} c_{\beta}^2 + B_{1} c_{\beta} + B_{0} \\
f_{Landau}(\eta_1, \eta_2, \eta_3) = (\eta_1^2 + \eta_2^2 + \eta_3^2) - 2(\eta_1^3 + \eta_2^3 + \eta_3^3) +  (\eta_1^4 + \eta_2^4 + \eta_3^4) + 5 (\eta_1^2 \eta_2^2 + \eta_2^2 \eta_3^2 + \eta_1^2 \eta_3^2) +  5(\eta_1^2 \eta_2^2 \eta_3^2) \\
H(\eta_p) = 3 \eta_p^2 - 2 \eta_p^3
\end{gather}
