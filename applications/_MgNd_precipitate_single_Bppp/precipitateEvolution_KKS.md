# KKS Phase Field Model of Precipitate Evolution (September 14th, 2023)
 
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
f_{grad}(\eta_1, \eta_2, \eta_3) = \frac{1}{2} \sum_{p=1}^3 \kappa_{ij}^{\eta_p} \eta_{p,i}  \eta_{p,j} 
\end{equation}
$$

$$
\begin{gather}
f_{elastic}(c,\eta_1, \eta_2, \eta_3,\epsilon) = \frac{1}{2} C_{ijkl}(\eta_1, \eta_2, \eta_3)  \left( \varepsilon_{ij} - \varepsilon ^0_{ij}(c, \eta_1, \eta_2, \eta_3) \right)\left( \varepsilon_{kl} - \varepsilon^0_{kl}(c, \eta_1, \eta_2, \eta_3)\right) 
\end{gather}
$$

$$
\begin{gather}
\varepsilon^0(c, \eta_1, \eta_2, \eta_3) = H(\eta_1) \varepsilon^0_{\eta_1} (c_{\beta})+ H(\eta_2) \varepsilon^0_{\eta_2} (c_{\beta}) + H(\eta_3) \varepsilon^0_{\eta_3} (c_{\beta}) 
C(\eta_1, \eta_2, \eta_3) = H(\eta_1) C_{\eta_1}+ H(\eta_2) C_{\eta_2} + H(\eta_3) C_{\eta_3} + \left( 1- H(\eta_1)-H(\eta_2)-H(\eta_3)\right)  C_{\alpha}
\end{gather}
$$

Here $\varepsilon^0_{\eta_p}$ are the composition dependent stress free strain transformation tensor corresponding to each structural order parameter, which is a function of the $\beta$ phase concentration, $c_{\beta}$, defined below.

In the KKS model (Kim 1999), the interfacial region is modeled as a mixture of the $\alpha$ and $\beta$ phases with concentrations $c_{alpha}$ and $c_{beta}$, respectively. The homogenous free energies for each phase, $f_{\alpha}$ and $f_{\beta}$ in this case, are typically given as functions of $c_{\alpha}$ and $c_{\beta}$, rather than directly as functions of $c$ and $\eta_p$. Thus, $f_{chem}(c, \eta_1, \eta_2, \eta_3)$ can be rewritten as 

$$
\begin{equation}
f_{chem}(c, \eta_1, \eta_2, \eta_3) = f_{\alpha}(c_{\alpha}) \left( 1- \sum_{p=1}^3 H(\eta_p)\right) + f_{\beta}(c_{\beta}) \sum_{p=1}^3 H(\eta_p)+ W f_{Landau}(\eta_1, \eta_2, \eta_3) 
\end{equation}
$$

The concentration in each phase is determined by the following system of equations:

$$
\begin{gather}
c =  c_{\alpha} \left( 1- \sum_{p=1}^3 H(\eta_p)\right) + c_{\beta} \sum_{p=1}^3 H(\eta_p) \\
\frac{\partial f_{\alpha}(c_{\alpha})}{\partial c_{\alpha}} = \frac{\partial f_{\beta}(c_{\beta})}{\partial c_{\beta}}
\end{gather}
$$

Given the following parabolic functions for the single-phase homogenous free energies:

$$
\begin{gather}
f_{\alpha}(c_{\alpha}) = A_{2} c_{\alpha}^2 + A_{1} c_{\alpha} + A_{0} \\
f_{\beta}(c_{\beta}) = B_{2} c_{\beta}^2 + B_{1} c_{\beta} + B_{0}
\end{gather}
$$

the single-phase concentrations are:

$$
\begin{gather}
c_{\alpha} = \frac{ B_2 c + \frac{1}{2} (B_1 - A_1) \sum_{p=1}^3 H(\eta_p) }{A_2 \sum_{p=1}^3 H(\eta_p) + B_2 \left( 1- \sum_{p=1}^3 H(\eta_p)\right) } \\
c_{\beta} =  \frac{ A_2 c + \frac{1}{2} (A_1 - B_1) \left[1-\sum_{p=1}^3 H(\eta_p)\right] }{A_2 \sum_{p=1}^3 H(\eta_p) + B_2 \left[ 1- \sum_{p=1}^3 H(\eta_p)\right] } 
\end{gather}
$$

## Required inputs

- $f_{\alpha}(c_{\alpha}), f_{\beta}(c_{\beta})$ - Homogeneous chemical free energy of the components of the binary system, example form given above
- $f_{Landau}(\eta_1, \eta_2, \eta_3)$ - Landau free energy term that controls the interfacial energy and prevents precipitates with different orientation varients from overlapping, example form given in Appendix I
- \$W$ - Barrier height for the Landau free energy term, used to control the thickness of the interface 
- $H(\eta_p)$ - Interpolation function for connecting the $\alpha$ phase and the $p^{th}$ orientation variant of the $\beta$ phase, example form given in Appendix I
- $\kappa^{\eta_p}$  - gradient penalty tensor for the $p^{th}$ orientation variant of the $\beta$ phase
- $C_{\eta_p}$ - fourth order elasticity tensor (or its equivalent second order Voigt representation) for the $p^{th}$ orientation variant of the $\beta$ phase
- $C_{\alpha}$ - fourth order elasticity tensor (or its equivalent second order Voigt representation) for the $\alpha$ phase
- $\varepsilon^0_{\eta_p}$ - stress free strain transformation tensor for the $p^{th}$ orientation variant of the $\beta$ phase


In addition, to drive the kinetics, we need:
- $M$  - mobility value for the concentration field
- $L$  - mobility value for the structural order parameter field


## Variational treatment
We obtain chemical potentials for the chemical potentials for the concentration and the structural order parameters by taking variational derivatives of $\Pi$:

$$
\begin{align}
  \mu_{c}  &= f_{\alpha,c} \left( 1- H(\eta_1)-H(\eta_2)-H(\eta_3)\right) +f_{\beta,c} \left(  H(\eta_1)  + H(\eta_2) + H(\eta_3) \right)  + C_{ijkl} (- \varepsilon^0_{ij,c}) \left( \varepsilon_{kl} - \varepsilon^0_{kl}\right) 
\end{align}
$$

$$
\begin{align}
\mu_{\eta_p}  &= [ f_{\beta}-f_{\alpha} -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}}] H(\eta_p){,\eta_p} + W f_{Landau,\eta_p}-
\kappa_{ij}^{\eta_p} \eta_{p,ij} + C_{ijkl} (- \varepsilon^0_{ij,\eta_p}) \left( \varepsilon_{kl} - \varepsilon^0_{kl}\right) + \frac{1}{2} C_{ijkl,\eta_p} \left( \varepsilon_{ij} - \varepsilon_{ij}^0 \right) \left( \varepsilon_{kl} - \varepsilon_{kl}^0\right)
\end{align}
$$

## Kinetics
Now the PDE for Cahn-Hilliard dynamics is given by:

$$
\begin{align}
\frac{\partial c}{\partial t} &= ~\nabla \cdot \left( \frac{1}{f_{,cc}}M \nabla \mu_c \right)
\end{align}
$$
 
where $M$ is a constant mobility and the factor of $\frac{1}{f_{,cc}}$ is added to guarentee constant diffusivity in the two phases. The PDE for Allen-Cahn dynamics is given by:

$$
\begin{align}
\frac{\partial \eta_p}{\partial t} &= - L \mu_{\eta_p} 
\end{align}
$$

where  $L$ is a constant mobility. 

## Mechanics
Considering variations on the displacement $u$ of the from $u+\epsilon w$, we have

$$
\begin{align}
\delta_u \Pi &=  \int_{\Omega}   \nabla w :  C(\eta_1, \eta_2, \eta_3) : \left( \varepsilon - \varepsilon^0(c,\eta_1, \eta_2, \eta_3)\right) ~dV = 0 
\end{align}
$$

where $\sigma = C(\eta_1, \eta_2, \eta_3) : \left( \varepsilon - \varepsilon^0(c,\eta_1, \eta_2, \eta_3)\right)$ is the stress tensor.

## Time discretization
Using forward Euler explicit time stepping, equations 

$$
\begin{align}
\frac{\partial c}{\partial t} &= ~\nabla \cdot \left( \frac{1}{f_{,cc}}M \nabla \mu_c \right)
\end{align}
$$

and 

$$
\begin{align}
\frac{\partial \eta_p}{\partial t} &= - L \mu_{\eta_p} 
\end{align}
$$

become:

$$
\begin{align}
c^{n+1} = c^{n}+\Delta t \left[\nabla \cdot \left(\frac{1}{f_{,cc}} M \nabla \mu_c \right) \right]
\end{align}
$$

$$
\begin{align}
\eta_p^{n+1} = \eta_p^n -\Delta t L \mu_{\eta_p}
\end{align}
$$

## Weak formulation
Writing equations 

$$
\begin{align}
\frac{\partial c}{\partial t} &= ~\nabla \cdot \left( \frac{1}{f_{,cc}}M \nabla \mu_c \right)
\end{align}
$$

and 

$$
\begin{align}
\frac{\partial \eta_p}{\partial t} &= - L \mu_{\eta_p} 
\end{align}
$$

in the weak form, with the arbirary variation given by $w$ yields:

$$
\begin{align}
\int_\Omega w c^{n+1} dV &= \int_\Omega wc^{n}+w  \Delta t \left[\nabla \cdot \left(\frac{1}{f_{,cc}}  M \nabla \mu_c \right) \right] dV
\end{align}
$$

$$
\begin{align}
%&= \int_\Omega wc^{n}+\nabla w \cdot (\Delta t  M \nabla \mu_c ) dV 
\end{align}
$$

$$
\begin{align}
\int_\Omega w \eta_p^{n+1} dV &= \int_\Omega w \eta_p^{n}-w  \Delta t L \mu_{\eta_p} dV
\end{align}
$$

The expression of $\frac{1}{f_{,cc}} \mu_c$ can be written as:

$$
\begin{align}
\frac{1}{f_{,cc}}  \nabla \mu_c = & \nabla c + (c_{\alpha}-c_{\beta}) \sum_{p=1}^3 H(\eta_p)_{,\eta_p} \nabla \eta_p 
\end{align}
$$

$$
\begin{align}
&+ \frac{1}{f_{,cc}} \left(\sum_{p=1}^3 (C_{ijkl}^{\eta_p} - C_{ijkl}^{\alpha} )\nabla \eta_p H_{,\eta_p}(\eta_p) \right)(-\epsilon_{ij,c}^0)(\epsilon_{ij} - \epsilon_{ij}^0)
\end{align}
$$

$$
\begin{align}
&- \frac{1}{f_{,cc}} C_{ijkl} \left(  \sum_{p=1}^3 \left( H_{,\eta_p}(\eta_p) \epsilon_{ij,c}^{0\eta_p} + \sum_{q=1}^3 \left( H(\eta_p) \epsilon_{ij,c\eta_q}^{0\eta_p} \right) \right) \nabla \eta_p + H(\eta_p) \epsilon_{ij,cc}^{0\eta_p} \nabla c \right)(\epsilon_{kl}-\epsilon_{kl}^0)
\end{align}
$$

$$
\begin{align}
&+ \frac{1}{f_{,cc}} C_{ijkl} (-\epsilon_{ij,c}^0) \left( \nabla \epsilon_{kl} -  \left( \sum_{p=1}^3 \left(H_{,\eta_p}(\eta_p) \epsilon_{kl}^{0\eta_p} -\sum_{q=1}^3 \epsilon_{kl,\eta_q}^{\eta_q} H(\eta_q) \right)\nabla \eta_p + H(\eta_p) \epsilon_{kl,c}^{0\eta_p} \nabla c \right) \right)
\end{align}
$$

Applying the divergence theorem to equation

$$
\begin{align}
\int_\Omega w c^{n+1} dV &= \int_\Omega wc^{n}+w  \Delta t \left[\nabla \cdot \left(\frac{1}{f_{,cc}}  M \nabla \mu_c \right) \right] dV
\end{align}
$$

one can derive the residual terms $r_c$ and $r_{cx}$:

$$
\begin{equation}
\int_\Omega w c^{n+1} dV = \int_\Omega wc^{n} +\nabla w \cdot (-\Delta t  M \frac{1}{f_{,cc}} \nabla \mu_c) dV
\end{equation}
$$

$$
\begin{align}
r_c &= c^{n}
\end{align}
$$

$$
\begin{align}
r_{cx} &= -\Delta t  M \frac{1}{f_{,cc}} \nabla \mu_c
\end{align}
$$

Expanding $\mu_{\eta_p}$ in equation 

$$
\begin{align}
\int_\Omega w \eta_p^{n+1} dV &= \int_\Omega w \eta_p^{n}-w  \Delta t L \mu_{\eta_p} dV
\end{align}
$$

and applying the divergence theorem yields the residual terms $r_{\eta_p}$ and $r_{\eta_p x}$:

$$
\begin{align}
\int_\Omega w \eta_p^{n+1} dV &=
\end{align}
$$

$$
\begin{align}
&\int_\Omega w \left(\eta_p^{n}-\Delta t L \left( (f_{\beta}-f_{\alpha})H_{,\eta_p}(\eta_p^n) -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}}H_{,\eta_p}(\eta_p^n) + W f_{Landau,\eta_p}
-C_{ijkl}  (H_{,\eta_p}(\eta_p) \epsilon_{ij}^{0 \eta_p}) (\epsilon_{kl} - \epsilon_{kl}^{0}) + \frac{1}{2} \left((C_{ijkl}^{\eta_p} - C_{ijkl}^{\alpha}) H_{,\eta_p}(\eta_p) \right) (\epsilon_{ij} - \epsilon_{ij}^{0}) (\epsilon_{kl} - \epsilon_{kl}^{0}) \right) \right) &+ \nabla w \cdot \left(-\Delta t  L \kappa_{ij}^{\eta_p} \eta_{p,i}^n \right) dV 
\end{align}
$$

where 

$$
\begin{align}
r_{\eta_p} &= \eta_p^{n}-\Delta t L \left( (f_{\beta}-f_{\alpha})H_{,\eta_p}(\eta_p^n) -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}}H_{,\eta_p}(\eta_p^n) + W f_{Landau,\eta_p}
-C_{ijkl}  (H_{,\eta_p}(\eta_p) \epsilon_{ij}^{0 \eta_p}) (\epsilon_{kl} - \epsilon_{kl}^{0}) + \frac{1}{2} \left((C_{ijkl}^{\eta_p} - C_{ijkl}^{\alpha}) H_{,\eta_p}(\eta_p) \right) (\epsilon_{ij} - \epsilon_{ij}^{0}) (\epsilon_{kl} - \epsilon_{kl}^{0}) \right)
\end{align}
$$

$$
\begin{align}
r_{\eta_p x} &= -\Delta t  L \kappa_{ij}^{\eta_p} \eta_{p,i}^n
\end{align}
$$

## Appendix I: Example functions for $f_{\alpha}$, $f_{\beta}$, $f_{Landau}$, $H(\eta_p)$ 

$$
\begin{gather}
f_{\alpha}(c_{\alpha}) = A_{2} c_{\alpha}^2 + A_{1} c_{\alpha} + A_{0} \\
f_{\beta}(c_{\beta}) = B_{2} c_{\beta}^2 + B_{1} c_{\beta} + B_{0} \\
f_{Landau}(\eta_1, \eta_2, \eta_3) = (\eta_1^2 + \eta_2^2 + \eta_3^2) - 2(\eta_1^3 + \eta_2^3 + \eta_3^3) +  (\eta_1^4 + \eta_2^4 + \eta_3^4) + 5 (\eta_1^2 \eta_2^2 + \eta_2^2 \eta_3^2 + \eta_1^2 \eta_3^2) +  5(\eta_1^2 \eta_2^2 \eta_3^2) \\
H(\eta_p) = 3 \eta_p^2 - 2 \eta_p^3
\end{gather}
$$
