
# PRISMS-PF Application Formulation:
 
\section{Variational formulation}
The total free energy of the system (neglecting boundary terms) is of the form,
\begin{equation}
\Pi(c, \eta_1, \eta_2, \eta_3, \Bepsilon) = \int_{\Omega} f(c, \eta_1, \eta_2, \eta_3, \Bepsilon) ~dV 
\end{equation}
where $c$ is the concentration of the $\beta$ phase, $\eta_p$ are the structural order parameters and $\Bvarepsilon$ is the small strain tensor. $f$, the free energy density is given by
\begin{equation}
 f(c, \eta_1, \eta_2, \eta_3, \Bepsilon) =   f_{chem}(c, \eta_1, \eta_2, \eta_3) + f_{grad}(\eta_1, \eta_2, \eta_3) + f_{elastic}(c,\eta_1, \eta_2, \eta_3,\Bepsilon)
\end{equation}
where
\begin{gather}
f_{chem}(c, \eta_1, \eta_2, \eta_3) = f_{\alpha}(c) \left( 1- H(\eta_1)-H(\eta_2)-H(\eta_3)\right) + f_{\beta}(c) \left(  H(\eta_1)  + H(\eta_2) + H(\eta_3)) \right) \\
f_{grad}(\eta_1, \eta_2, \eta_3) = \frac{1}{2} \sum_{p=1}^3 \Bkappa^{\eta_p}_{ij} \eta_{p,i}  \eta_{p,j} \\
f_{elastic}(c,\eta_1, \eta_2, \eta_3,\Bepsilon) = \frac{1}{2} \bC_{ijkl}(\eta_1, \eta_2, \eta_3)  \left( \Bvarepsilon_{ij} - \Bvarepsilon ^0_{ij}(c, \eta_1, \eta_2, \eta_3) \right)\left( \Bvarepsilon_{kl} - \Bvarepsilon^0_{kl}(c, \eta_1, \eta_2, \eta_3)\right) \\
\Bvarepsilon^0(c, \eta_1, \eta_2, \eta_3) = H(\eta_1) \Bvarepsilon^0_{\eta_1} (c)+ H(\eta_2) \Bvarepsilon^0_{\eta_2} (c) + H(\eta_3) \Bvarepsilon^0_{\eta_3} (c) \\
\bC(\eta_1, \eta_2, \eta_3) = H(\eta_1) \bC_{\eta_1}+ H(\eta_2) \bC_{\eta_2} + H(\eta_3) \bC_{\eta_3} + \left( 1- H(\eta_1)-H(\eta_2)-H(\eta_3)\right)  \bC_{\alpha}
\end{gather}
Here $\Bvarepsilon^0_{\eta_p}$ are the composition dependent stress free strain transformation tensor corresponding to each structural order parameter.

\section{Required inputs}
\begin{itemize}
\item $f_{\alpha}(c), f_{\beta}(c)$ - Homogeneous chemical free energy of the components of the binary system, example form given in Appendix I
\item $H(\eta_p)$ - Interpolation function for connecting the $\alpha$ phase and the $p^{th}$ orientation variant of the $\beta$ phase, example form given in Appendix I
\item $\Bkappa^{\eta_p}$  - gradient penalty tensor for the $p^{th}$ orientation variant of the $\beta$ phase
\item $\bC_{\eta_p}$ - fourth order elasticity tensor (or its equivalent second order Voigt representation) for the $p^{th}$ orientation variant of the $\beta$ phase
\item $\bC_{\alpha}$ - fourth order elasticity tensor (or its equivalent second order Voigt representation) for the $\alpha$ phase
\item $\Bvarepsilon^0_{\eta_p}$ - stress free strain transformation tensor for the $p^{th}$ orientation variant of the $\beta$ phase
\end{itemize}
In addition, to drive the kinetics, we need:
\begin{itemize}
\item $M$  - mobility value for the concentration field
\item $L$  - mobility value for the structural order parameter field
\end{itemize}

\section{Variational treatment}
From the variational derivatives given in Appendix II, we obtain the chemical potentials for the concentration and the structural order parameters:
\begin{align}
  \mu_{c}  &= f_{\alpha,c} \left( 1- H(\eta_1)-H(\eta_2)-H(\eta_3)\right) +f_{\beta,c} \left(  H(\eta_1)  + H(\eta_2) + H(\eta_3) \right)  + \bC_{ijkl} (- \Bvarepsilon^0_{ij,c}) \left( \Bvarepsilon_{kl} - \Bvarepsilon^0_{kl}\right) \\
  \mu_{\eta_p}  &= (f_{\beta}-f_{\alpha})H(\eta_p)_{,\eta_p} - \Bkappa^{\eta_p}_{ij} \eta_{p,ij} + \bC_{ijkl} (- \Bvarepsilon^0_{ij,\eta_p}) \left( \Bvarepsilon_{kl} - \Bvarepsilon^0_{kl}\right) + \frac{1}{2} \bC_{ijkl,\eta_p} \left( \Bvarepsilon_{ij} - \Bvarepsilon ^0_{ij} \right) \left( \Bvarepsilon_{kl} - \Bvarepsilon^0_{kl}\right)
\end{align}

\section{Kinetics}
Now the PDE for Cahn-Hilliard dynamics is given by:
\begin{align}
  \frac{\partial c}{\partial t} &= ~\grad \cdot (M \grad \mu_c) \label{CH_eqn}
  \end{align}
  and the PDE for Allen-Cahn dynamics is given by:
  \begin{align}
    \frac{\partial \eta_p}{\partial t} &= - L \mu_{\eta_p} \label{AC_eqn}
\end{align}
where $M$ and $L$ are the constant mobilities. 

\section{Mechanics}
Considering variations on the displacement $u$ of the from $u+\epsilon w$, we have
\begin{align}
\delta_u \Pi &=  \int_{\Omega}   \grad w :  \bC(\eta_1, \eta_2, \eta_3) : \left( \Bvarepsilon - \Bvarepsilon^0(c,\eta_1, \eta_2, \eta_3)\right) ~dV = 0 \\
\end{align}
where $\Bsigma = \bC(\eta_1, \eta_2, \eta_3) : \left( \Bvarepsilon - \Bvarepsilon^0(c,\eta_1, \eta_2, \eta_3)\right)$ is the stress tensor. \\

Now consider\\
\begin{align}
R &=  \int_{\Omega}   \grad w :  \bC(\eta_1, \eta_2, \eta_3) : \left( \Bvarepsilon - \Bvarepsilon^0(c,\eta_1, \eta_2, \eta_3)\right) ~dV = 0 
\end{align}
We solve for $R=0$ using a gradient scheme which involves the following linearization:
\begin{align}
R~|_{u}+ \frac{\delta R}{\delta u} \Delta u &= 0 \\
\Rightarrow \frac{\delta R}{\delta u} \Delta u &= -R~|_{u}
\end{align}
This is the linear system $Ax=b$ which we solve implicitly using the Conjugate Gradient scheme. For clarity, here in the left hand side (LHS) $A=\frac{\delta R}{\delta u}$, $x=\Delta u$ and the right hand side (RHS) is $b=-R~|_{u}$.


\section{Time discretization}
Using forward Euler explicit time stepping, equations \ref{CH_eqn} and \ref{AC_eqn} become:
\begin{align}
c^{n+1} = c^{n}+\Delta t [\nabla \cdot (M \nabla \mu_c) ]\\
\eta_p^{n+1} = \eta_p^n -\Delta t L \mu_{\eta_p}
\end{align}

\section{Weak formulation and residual expressions}
\subsection{The Cahn-Hillard and Allen-Cahn equations}
Writing equations \ref{CH_eqn} and \ref{AC_eqn} in the weak form, with the arbirary variation given by $w$ yields:
\begin{align}
\int_\Omega w c^{n+1} dV &= \int_\Omega wc^{n}+w  \Delta t [\nabla \cdot (M \nabla \mu_c) ] dV \label{CH_weak} \\
%&= \int_\Omega w\underbrace{c^{n}}_{r_c}+\nabla w \cdot (\underbrace{\Delta t  M \nabla \mu_c}_{r_{cx}} ) dV \\
\int_\Omega w \eta_p^{n+1} dV &= \int_\Omega w \eta_p^{n}-w  \Delta t L \mu_{\eta_p} dV  \label{AC_weak}
%&= \int_\Omega w\underbrace{c^{n}}_{r_c}+\nabla w \cdot (\underbrace{\Delta t  M \nabla \mu_c}_{r_{cx}} ) dV 
\end{align}

The gradient of $\mu_c$ is:
\begin{equation}
\begin{split}
\nabla \mu_c = & \nabla c \left[f_{\alpha,cc}+ \sum_{p=1}^3 H(\eta_p)(f_{\beta,cc}-f_{\alpha,cc}) \right]+ \sum_{p=1}^3 \nabla \eta_p H(\eta_p)_{,\eta_p} (f_{\beta,c}-f_{\alpha,c}) \\
&+ \left[ \sum_{p=1}^3 (C_{ijkl}^{\eta_p} - C_{ijkl}^{\alpha} )\nabla \eta_p H(\eta_p)_{,\eta_p} \right](-\epsilon_{ij,c}^0)(\epsilon_{ij} - \epsilon_{ij}^0) \\
&- C_{ijkl} \left[  \sum_{p=1}^3 H(\eta_p)_{,\eta_p} \epsilon_{ij,c}^{0\eta_p} \nabla \eta_p + H(\eta_p) \epsilon_{ij,cc}^{0\eta_p} \nabla c \right](\epsilon_{kl}-\epsilon_{kl}^0)\\
&+ C_{ijkl} (-\epsilon_{ij,c}^0) \left[ \nabla \epsilon_{ij} -  \left( \sum_{p=1}^3 H(\eta_p)_{,\eta_p} \epsilon_{kl}^{0\eta_p} \nabla \eta_p + H(\eta_p) \epsilon_{kl,c}^{0\eta_p} \nabla c \right) \right]
\end{split}
\end{equation}

Applying the divergence theorem to equation \ref{CH_weak}, one can derive the residual terms $r_c$ and $r_{cx}$:
\begin{equation}
\int_\Omega w c^{n+1} dV = \int_\Omega w\underbrace{c^{n}}_{r_c}+\nabla w \cdot (\underbrace{-\Delta t  M \nabla \mu_c}_{r_{cx}} ) dV
\end{equation}

Expanding $\mu_{\eta_p}$ in equation \ref{AC_weak} and applying the divergence theorem yields the residual terms $r_{\eta_p}$ and $r_{\eta_p x}$:
\begin{equation}
\begin{split}
\int_\Omega w \eta_p^{n+1} dV &= \\
&\int_\Omega w \Bigg\{\underbrace{\eta_p^{n}-\Delta t L \bigg[(f_{\beta}-f_{\alpha})H(\eta_p^n)_{,\eta_p} - C_{ijkl} \left( H(\eta_p)_{,\eta_p} \epsilon_{ij}^{0 \eta_p}\right)\left(\epsilon_{kl} - \epsilon_{kl}^{0} \right) }_{r_{\eta_p}} \\
&\underbrace{+ \frac{1}{2} \left[ (C_{ijkl}^{\eta_p} - C_{ijkl}^{\alpha}) H(\eta_p)_{,\eta_p} \right] \left(\epsilon_{ij} - \epsilon_{ij}^{0} \right) \left(\epsilon_{kl} - \epsilon_{kl}^{0} \right) \bigg] \Bigg\}}_{r_{\eta_p}~cont.}  \\
&+ \nabla w \cdot (\underbrace{-\Delta t  L \Bkappa^{\eta_p}_{ij} \eta_{p,i}^n}_{r_{\eta_p x}} ) dV 
\end{split}
\end{equation}

The above values of $r_c$, $r_{cx}$,  $r_{\eta_p}$, and $r_{\eta_p x}$are used to define the residuals in the following input file: \\
$applications/precipitateEvolution/equations.h$

\subsection{The mechanical equilbrium equation}
In PRISMS-PF, two sets of residuals are required for elliptic PDEs (such as this one), one for the left-hand side of the equation (LHS) and one for the right-hand side of the equation (RHS). We solve $R=\delta_u \Pi$ by casting this in a form that can be solved as a matrix inversion problem. This will involve a brief detour into the discretized form of the equation. First we derive an expression for the solution, given an initial guess, $u_0$:
\begin{gather}
0 = R(u) = R(u_0 + \Delta u)
\end{gather}
where $\Delta u = u - u_0$. Then, applying the discretization that $u = \sum_i w^i U^i$, we can write the following linearization:
\begin{equation}
\frac{\delta R(u)}{\delta u} \Delta U = -R(u_0) \label{matrix_eqn}
\end{equation}
The discretized form of this equation can be written as a matrix inversion problem. However, in PRISMS-PF, we only care about the product $\frac{\delta R(u)}{\delta u} \Delta U$. Taking the variational derivative of $R(u)$ yields:
\begin{align}
\frac{\delta R(u)}{\delta u} &= \frac{d}{d\alpha} \int_{\Omega}   \nabla w :C: \left[ \epsilon (u+\alpha w) - \epsilon^0 \right] ~dV  \bigg{|}_{\alpha=0} \\
&=  \int_{\Omega}   \nabla w :C: \frac{1}{2}\frac{d}{d\alpha}\left[ \nabla(u+\alpha w) + \nabla(u+\alpha w)^T  - \epsilon^0\right] ~dV \bigg{|}_{\alpha=0}\\
&= \int_{\Omega}   \nabla w :C: \frac{d}{d\alpha} \left[ \nabla(u+\alpha w) - \epsilon^0 \right]  ~dV \bigg{|}_{\alpha=0} \quad (due ~to ~the ~symmetry ~of ~C) \\
&= \int_{\Omega}   \nabla w :C: \nabla w  ~dV 
\end{align}
In its discretized form $\frac{\delta R(u)}{\delta u} \Delta U$ is:
\begin{equation}
\frac{\delta R(u)}{\delta u} \Delta U = \sum_i \sum_j \int_{\Omega} \nabla N^i : C : \nabla N^j dV ~\Delta U^j
\end{equation}
Moving back to the non-discretized form yields:
\begin{equation}
\frac{\delta R(u)}{\delta u} \Delta U = \int_{\Omega} \nabla w : C : \nabla (\Delta u) dV
\end{equation}
Thus, the full equation relating $u_0$ and $\Delta u$ is:
\begin{equation}
\int_{\Omega} \nabla w : \underbrace{C : \nabla (\Delta u)}_{r_{ux}^{LHS}} dV = -\int_{\Omega}   \grad w : \underbrace{\sigma}_{r_{ux}} ~dV
\end{equation}
The above values of $r_{ux}^{LHS}$ and $r_{ux}$ are used to define the residuals in the following input file: \\
$applications/precipitateEvolution/equations.h$

\section{Appendix I: Example functions for $f_{\alpha}$, $f_{\beta}$, $H(\eta_p)$ }
\begin{gather}
f_{\alpha}(c) = A_{2, \alpha} c^2 + A_{1, \alpha} c + A_{0, \alpha} \\
f_{\beta}(c) = A_{2, \beta} c^2 + A_{1, \beta} c + A_{0, \beta} \\
H(\eta_p) = 10 \eta_p^3 - 15 \eta_p^4 + 6 \eta_p^5
\end{gather}

\section{Appendix II: Variational Derivatives}
Variational derivative of $\Pi$ with respect to $\eta_p$ (where $\eta_q$ and $\eta_r$ correspond to the structural order parameters for the other two orientational variants):
\begin{gather}
\delta_{\eta_p} \Pi  =  \frac{d}{d\alpha} \left[\int_{\Omega}  f_{chem}(c,\eta_p+\alpha w,\eta_q,\eta_r) + f_{grad}(\eta_p+\alpha w,\eta_q,\eta_r) + f_{el}(c,\eta_p+\alpha w,\eta_q,\eta_r,\epsilon) dV  \right]_{\alpha=0}
\end{gather}
Breaking up each of these terms yields:

\begin{align}
\begin{split}
\frac{d}{d\alpha} \left[ f_{chem}(c,\eta_p+\alpha w,\eta_q,\eta_r)\right]_{\alpha=0} &= f_{\alpha}(c) \left[  -\frac{\partial H(\eta_p+\alpha w)}{\partial (\eta_p + \alpha w)} \frac{\partial(\eta_p + \alpha w)}{\partial \alpha} \right]_{\alpha=0} \\
&+f_{\beta}(c)  \left[  \frac{\partial H(\eta_p+\alpha w)}{\partial (\eta_p + \alpha w)} \frac{\partial(\eta_p + \alpha w)}{\partial \alpha} \right]_{\alpha=0} \\
\\
&=f_{\alpha}(c) \left[  -\frac{\partial H(\eta_p)}{\partial \eta_p} w \right] 
+f_{\beta}(c) \left[  \frac{\partial H(\eta_p)}{\partial \eta_p } w \right] 
\end{split}
\end{align}

\begin{align}
\begin{split}
\frac{d}{d\alpha} \left[ f_{grad}(\eta_p+\alpha w,\eta_q,\eta_r)\right]_{\alpha=0} &= \frac{1}{2} \left[ \kappa_{ij}^{\eta_p} (\eta_p+\alpha w)_{,i}(\eta_p+\alpha w)_{,j} +\kappa_{ij}^{\eta_q} (\eta_q)_{,i}(\eta_q)_{,j} + \kappa_{ij}^{\eta_r} (\eta_r)_{,i}(\eta_r)_{,j}  \right]_{\alpha=0} \\
\\
&= \kappa_{ij} w_{,i} \eta_{p,j}
\end{split}
\end{align}

\begin{align}
\begin{split}
\frac{d}{d\alpha} \left[ f_{el}(c,\eta_p+\alpha w,\eta_q,\eta_r,\epsilon)\right]_{\alpha=0} &= \frac{1}{2}  \bigg[  \frac{\partial C_{ijkl}(\eta_p+\alpha w,\eta_q,\eta_r)}{\partial (\eta_p + \alpha w)} \frac{\partial(\eta_p + \alpha w)}{\partial \alpha} \\ 
&\cdot\big (\epsilon_{ij}-\epsilon_{ij}^0 (c,\eta_p+\alpha w,\eta_q,\eta_r)\big) \big(\epsilon_{kl}-\epsilon_{kl}^0 (c,\eta_p+\alpha w,\eta_q,\eta_r)\big) \\
&+ C_{ijkl}(\eta_p+\alpha w,\eta_q,\eta_r) \bigg (-\frac{\partial\epsilon_{ij}^0 (c,\eta_p+\alpha w,\eta_q,\eta_r)}{\partial (\eta_p + \alpha w)} \frac{\partial(\eta_p + \alpha w)}{\partial \alpha} \bigg) \\ &\cdot \big(\epsilon_{kl}-\epsilon_{kl}^0 (c,\eta_p+\alpha w,\eta_q,\eta_r)\big) \\
&+ C_{ijkl}(\eta_p+\alpha w,\eta_q,\eta_r)   \big(\epsilon_{ij}-\epsilon_{ij}^0 (c,\eta_p+\alpha w,\eta_q,\eta_r)\big) \\ 
&\cdot\bigg (-\frac{\partial\epsilon_{kl}^0 (c,\eta_p+\alpha w,\eta_q,\eta_r)}{\partial (\eta_p + \alpha w)} \frac{\partial(\eta_p + \alpha w)}{\partial \alpha} \bigg)\bigg]_{\alpha=0} \\
\\
&= \frac{1}{2}  \bigg[  \frac{\partial C_{ijkl}(\eta_p,\eta_q,\eta_r)}{\partial \eta_p} w \big(\epsilon_{ij}-\epsilon_{ij}^0 (c,\eta_p,\eta_q,\eta_r) \big) \big(\epsilon_{kl}-\epsilon_{kl}^0 (c,\eta_p,\eta_q,\eta_r) \big) \bigg] \\
&+ C_{ijkl}(\eta_p,\eta_q,\eta_r) \bigg (-\frac{\partial\epsilon_{ij}^0 (c,\eta_p,\eta_q,\eta_r)}{\partial \eta_p} w \bigg) \big(\epsilon_{kl}-\epsilon_{kl}^0 (c,\eta_p,\eta_q,\eta_r)\big) \\
\end{split}
\end{align}

Putting the terms back together yields:
\begin{gather}
\begin{split}
\delta_{\eta_p} \Pi  &=  \int_{\Omega}  f_{\alpha}(c) \left[  -\frac{\partial H(\eta_p)}{\partial \eta_p} w \right] +f_{\beta}(c) \left[  \frac{\partial H(\eta_p)}{\partial \eta_p } w \right] \\
&+ \kappa_{ij} w_{,i} \eta_{p,j} \\
&+\frac{1}{2}  \bigg[  \frac{\partial C_{ijkl}(\eta_p,\eta_q,\eta_r)}{\partial (\eta_p)} w \big(\epsilon_{ij}-\epsilon_{ij}^0 (c,\eta_p,\eta_q,\eta_r) \big) \big(\epsilon_{kl}-\epsilon_{kl}^0 (c,\eta_p,\eta_q,\eta_r) \big) \bigg] \\
&+ C_{ijkl}(\eta_p,\eta_q,\eta_r) \bigg (-\frac{\partial\epsilon_{ij}^0 (c,\eta_p,\eta_q,\eta_r)}{\partial \eta_p} w \bigg) \big(\epsilon_{kl}-\epsilon_{kl}^0 (c,\eta_p,\eta_q,\eta_r)\big)  ~dV  
\end{split}
\end{gather}

Variational derivative of $\Pi$ with respect to $c$ :
\begin{gather}
\delta_{c} \Pi  =  \frac{d}{d\alpha} \left[\int_{\Omega}  f_{chem}(c+\alpha w,\eta_p,\eta_q,\eta_r) + f_{grad}(\eta_p,\eta_q,\eta_r) + f_{el}(c+\alpha w,\eta_p,\eta_q,\eta_r,\epsilon) dV  \right]_{\alpha=0}
\end{gather}
Breaking up each of these terms yields:
\begin{align}
\begin{split}
\frac{d}{d\alpha} \left[ f_{chem}(c+\alpha w,\eta_p,\eta_q,\eta_r)\right]_{\alpha=0}  &= \bigg[ \frac{\partial f_{\alpha}(c + \alpha w)}{\partial(c+\alpha w)} \frac{\partial(c+\alpha w)}{\partial \alpha} \left(1-\sum_{p=1}^3 H(\eta_p)\right) \\
&+ \frac{\partial f_{\beta}(c+\alpha w)}{\partial(c+\alpha w)} \frac{\partial(c+\alpha w)}{\partial \alpha} \left(\sum_{p=1}^3 H(\eta_p)\right) \\
&= \frac{\partial f_{\alpha}(c)}{\partial c} w \left(1-\sum_{p=1}^3 H(\eta_p)\right) +\frac{\partial f_{\beta}(c)}{\partial c} w \left(\sum_{p=1}^3 H(\eta_p)\right)
\end{split}
\end{align}

\begin{align}
\frac{d}{d\alpha} \left[ f_{grad}(\eta_p,\eta_q,\eta_r)\right]_{\alpha=0} = 0
\end{align}

\begin{align}
\begin{split}
\frac{d}{d\alpha} \left[ f_{el}(c+\alpha w,\eta_p,\eta_q,\eta_r,\epsilon)\right]_{\alpha=0} &= \frac{1}{2} C_{ijkl}(\eta_p,\eta_q,\eta_r)  \bigg[ -\frac{\partial \epsilon_{ij}^0 (c+\alpha w,\eta_p,\eta_q,\eta_r)}{\partial (c + \alpha w)} \frac{\partial(c + \alpha w)}{\partial \alpha} \big(\epsilon_{kl}-\epsilon_{kl}^0 (c+\alpha w,\eta_p,\eta_q,\eta_r)\big)  \\ 
&- \big(\epsilon_{ij}-\epsilon_{ij}^0 (c+\alpha w,\eta_p,\eta_q,\eta_r)\big) \frac{\partial \epsilon_{ij}^0 (c+\alpha w,\eta_p,\eta_q,\eta_r)}{\partial (c + \alpha w)} \frac{\partial(c + \alpha w)}{\partial \alpha} \bigg]_{\alpha=0} \\
\\
&= -C_{ijkl}(\eta_p,\eta_q,\eta_r) \frac{\partial \epsilon_{ij}^0 (c,\eta_p,\eta_q,\eta_r)}{\partial c} w \big(\epsilon_{kl}-\epsilon_{kl}^0 (c+\alpha w,\eta_p,\eta_q,\eta_r)\big)
\end{split}
\end{align}

Putting the terms back together yields:
\begin{gather}
\begin{split}
\delta_{c} \Pi  &=  \int_{\Omega}  \frac{\partial f_{\alpha}(c)}{\partial c} w \left(1-\sum_{p=1}^3 H(\eta_p)\right) +\frac{\partial f_{\beta}(c)}{\partial c} w \left(\sum_{p=1}^3 H(\eta_p)\right) \\
&-C_{ijkl}(\eta_p,\eta_q,\eta_r) \frac{\partial \epsilon_{ij}^0 (c,\eta_p,\eta_q,\eta_r)}{\partial c} w \big(\epsilon_{kl}-\epsilon_{kl}^0 (c+\alpha w,\eta_p,\eta_q,\eta_r)\big) ~dV
\end{split}
\end{gather}
