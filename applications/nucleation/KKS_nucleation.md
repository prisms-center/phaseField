# KKS Phase Field Model of Precipitate Evolution coupled with nucleation (October 2, 2024)

The Nucleation Model application for PRISMS-PF  incorporates a stochastic method to add nuclei to the KKS phase field model for precipitate growth. Nuclei are seeded throughout the time evolution of the system based on a probability that depends on the local solute supersaturation. This document is divided in two sections. In the first section, the phase field model formulation for precipitate evolution in a binary alloy (without elastic effects) is presented. In the second section the nucleation method is presented.

## Precipitate Evolution
### Variational formulation
In the absence of elastic effects total free energy of the 2-component system (neglecting boundary terms) is of the form,

$$
\begin{equation}
\Pi(c, \eta) = \int_{\Omega} f(c, \eta) ~dV
\end{equation}
$$

where $c$ is the concentration of the $\beta$ phase and  $\eta$ is the set of structural order parameters. The free energy density, $f$, is given by

$$
\begin{equation}
 f(c, \eta) =   f_{chem}(c, \eta) + f_{grad}(\eta)
\end{equation}
$$

where

$$
\begin{equation}
f_{chem}(c, \eta) = f_{\alpha}(c,\eta) \left( 1- H(\eta)\right) + f_{\beta}(c,\eta) H(\eta)+ W f_{Landau}(\eta)
\end{equation}
$$

and

$$
\begin{equation}
f_{grad}(\eta) = \frac{1}{2} \kappa | \nabla \eta |^2 \\
\end{equation}
$$

In the KKS model (Kim 1999), the interfacial region is modeled as a mixture of the $\alpha$ and $\beta$ phases with concentrations $c_{\alpha}$ and $c_{\beta}$, respectively. The homogeneous free energies for each phase, $f_{\alpha}$ and $f_{\beta}$ in this case, are typically given as functions of $c_{\alpha}$ and $c_{\beta}$, rather than directly as functions of $c$ and $\eta_p$. Thus, $f_{chem}(c, \eta)$ can be rewritten as

$$
\begin{equation}
f_{chem}(c, \eta) = f_{\alpha}(c_\alpha) \left( 1- H(\eta)\right) + f_{\beta}(c_\beta) H(\eta)+ W f_{Landau}(\eta)
\end{equation}
$$

The concentration in each phase is determined by the following system of equations:

$$
\begin{align}
c =  c_{\alpha} \left( 1- H(\eta)\right) + c_{\beta} H(\eta)
\end{align}
$$

$$
\begin{align}
\frac{\partial f_{\alpha}(c_{\alpha})}{\partial c_{\alpha}} = \frac{\partial f_{\beta}(c_{\beta})}{\partial c_{\beta}}
\end{align}
$$

Given the following parabolic functions for the single-phase homogeneous free energies:

$$
\begin{align}
f_{\alpha}(c_{\alpha}) = A_{2} c_{\alpha}^2 + A_{1} c_{\alpha} + A_{0}
\end{align}
$$

$$
\begin{align}
f_{\beta}(c_{\beta}) = B_{2} c_{\beta}^2 + B_{1} c_{\beta} + B_{0}
\end{align}
$$

the single-phase concentrations are:

$$
\begin{align}
c_{\alpha} = \frac{ B_2 c + \frac{1}{2} (B_1 - A_1) H(\eta) }{A_2 H(\eta) + B_2 \left( 1- H(\eta)\right) }
\end{align}
$$

$$
\begin{align}
c_{\beta} =  \frac{ A_2 c + \frac{1}{2} (A_1 - B_1) \left[1-H(\eta)\right] }{A_2  H(\eta) + B_2 \left[ 1- H(\eta)\right] }
\end{align}
$$

### Required inputs

- $f_{\alpha}(c_{\alpha}), f_{\beta}(c_{\beta})$ - Homogeneous chemical free energy of the components of the binary system, example form given above
- $f_{Landau}(\eta)$ - Landau free energy term that controls the interfacial energy. Example form given in Appendix I
- $W$ - Barrier height for the Landau free energy term, used to control the thickness of the interface
- $H(\eta)$ - Interpolation function for connecting the $\alpha$ phase and the $\beta$ phase. Example form given in Appendix I
- $\kappa^{\eta_p}$  - gradient penalty coefficient for the $\alpha - \beta$ interface

In addition, to drive the kinetics, we need:

- item $M$  - mobility value for the concentration field
- item $L$  - mobility value for the structural order parameter field


### Variational treatment
We obtain chemical potentials for the concentration and the structural order parameter by taking variational derivatives of $\Pi$:

$$
\begin{align}
  \mu_{c}  &= f_{\alpha,c} \left( 1- H(\eta)\right) +f_{\beta,c} H(\eta)
\end{align}
$$

$$
\begin{align}
  \mu_{\eta}  &= \left[ f_{\beta}-f_{\alpha} -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}} \right] H_{,\eta}(\eta) + W f_{Landau,\eta}- \kappa\nabla^2\eta
\end{align}
$$

### Kinetics
Now the PDE for Cahn-Hilliard dynamics is given by:

$$
\begin{align}
  \frac{\partial c}{\partial t} &= ~\nabla \cdot \left( \frac{1}{f_{,cc}}M \nabla \mu_c \right)
  \end{align}
$$

  where $M$ is a constant mobility and the factor of $\frac{1}{f_{,cc}}$ is added to guarantee constant diffusivity in the two phases. The PDE for Allen-Cahn dynamics is given by:

$$
  \begin{align}
    \frac{\partial \eta}{\partial t} &= - L \mu_{\eta_p}
\end{align}
$$

where  $L$ is a constant mobility.

### Time discretization
Using forward Euler explicit time stepping, the equations from the Kinetics section become:

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

### Weak formulation
Writing the equations from the Kinetics section in the weak form, with the arbitrary variation given by $w$ yields:

$$
\begin{align}
\int_\Omega w c^{n+1} dV &= \int_\Omega wc^{n}+w  \Delta t \left[\nabla \cdot \left(\frac{1}{f_{,cc}}  M \nabla \mu_c \right) \right] dV
\end{align}
$$

$$
\begin{align}
%&= \int_\Omega wc^{n} +\nabla w \cdot (\Delta t  M \nabla \mu_c ) dV
\end{align}
$$

$$
\begin{align}
r_c &= c^{n}
\end{align}
$$

$$
\begin{align}
r_{cx} &= \Delta t  M \nabla \mu_c
\end{align}
$$

$$
\begin{align}
\int_\Omega w \eta^{n+1} dV &= \int_\Omega w \eta^{n}-w  \Delta t L \mu_{\eta} dV
\end{align}
$$

$$
\begin{align}
%&= \int_\Omega wc^{n} +\nabla w \cdot (\Delta t  M \nabla \mu_c ) dV
\end{align}
$$

$$
\begin{align}
r_c &= c^{n}
\end{align}
$$

$$
\begin{align}
r_{cx} &= \Delta t  M \nabla \mu_c
\end{align}
$$

The expression of $\frac{1}{f_{,cc}} \mu_c$ can be written as:

$$
\begin{equation}
\frac{1}{f_{,cc}}  \nabla \mu_c =  \nabla c + (c_{\alpha}-c_{\beta}) H(\eta)_{,\eta} \nabla \eta
\end{equation}
$$

Applying the divergence theorem to the weak CH equation, one can derive the residual terms $r_c$ and $r_{cx}$:

$$
\begin{equation}
\int_\Omega w c^{n+1} dV = \int_\Omega w c^{n} +\nabla w \cdot (-\Delta t  M \frac{1}{f_{,cc}} \nabla \mu_c ) dV
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

Expanding $\mu_{\eta}$ in the weak AC equation and applying the divergence theorem yields the residual terms $r_{\eta}$ and $r_{\eta x}$:

$$
\begin{align}
\int_\Omega w \eta^{n+1} dV = &\int_\Omega w \bigg[ \eta^{n}-\Delta t L \bigg[(f_{\beta}-f_{\alpha})H_{,\eta}(\eta^n) -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}}H_{,\eta}(\eta^n) + W f_{Landau,\eta}
\end{align}
$$

$$
\begin{align}
&+ \nabla w \cdot (-\Delta t  L \kappa \nabla \eta^n ) dV
\end{align}
$$

$$
\begin{align}
r_{\eta} &= \eta^{n}-\Delta t L \bigg[(f_{\beta}-f_{\alpha})H_{,\eta}(\eta^n) -(c_{\beta}-c_{\alpha}) f_{\beta,c_{\beta}}H_{,\eta}(\eta^n) + W f_{Landau,\eta}
\end{align}
$$

$$
\begin{align}
r_{\eta x} &= -\Delta t  L \kappa \nabla \eta^n
\end{align}
$$

## Nucleation method

We follow the same approach as Jokisaari and Thornton [Comput. Mater. Sci. **112**, 128-138 (2016)] which consists of adding nuclei throughout a phase field simulation based on a probability that depends on the local supersaturation. This probability is calculated every fixed number of time steps and for every element of the computational domain. In each nucleation event, nucleation is triggered at a point within the $\alpha$ phase. Each nucleus is then added to the system by modifying the order parameter to it's $\beta$ phase value within a small domain around the selected nucleation center. This domain can be spherical/circular or ellipsoidal/elliptical.

### Nucleation rate

From classical nucleation theory, the nucleation rate for critical nuclei $J^*$ is given by

$$
\begin{align}
J^* (\mathbf{r},t)=Zn\beta^* \exp \left( -\frac{\Delta G^* }{k_B T} \right) \exp \left( -\frac{\tau}{t} \right),
\end{align}
$$

where $Z$ is the Zeldovich factor, $n$ is the number of nucleation sites per volume, $\beta^*$ is the frequency at which a critical nucleus becomes supercritical, $\Delta G^*$ is the nucleation energy barrier, $k_B$ is the Boltzmann constant, $T$ is the temperature, $t$ is time and $\tau$ is the incubation time. It can be shown that, in the dilute limit and for constant temperature, the previous equation can be simplified by grouping approximately constant terms in both the exponential and pre-exponential factors:

$$
\begin{align}
J^*(\mathbf{r},t)=k_1\exp \left( -\frac{k_2}{(\Delta c)^{d-1}} \right) \exp \left(-\frac{\tau}{t} \right),
\end{align}
$$

where  $k_1$ and $k_2$ are now taken as constant parameters, $\Delta c=c(\mathbf{r},t)-c_\alpha^{eq}$ is the local supersaturation in the $\alpha$ phase and $d$ is the dimensionality of the system (*e.g.* $d=2$ or $d=3$).

### Nucleation probability

Considering  $J^*$ to be approximately constant within a small volume, $\Delta V$, and for a small time interval, $\Delta t$, the probability that at least one nucleation event occurs in $\Delta V$ within $\Delta t$ is given by [Simmons et al., Scripta Mater. **43**, 935 (2000)]

$$
\begin{align}
P(\mathbf{r},t) = 1 - \exp \left( -J^* \Delta V \Delta t \right)
\end{align}
$$

### Hold time

After each nucleus is added, there is a `hold' time interval, $\Delta t_h$, during which the order parameter value is fixed within a small window that encompasses the new nucleus. The purpose of this hold time is to allow the concentration to evolve within the nucleus to a value close to the coexistence composition for $\beta$ phase, and therefore, to create small a solute depleted zone around the nucleus. After the hold time, the nucleus is allowed to evolve into a precipitate.

### Required nucleation inputs

- $k_1$ - Constant pre-exponential factor in the equation for simple nucleation rate
- $k_2$ - Parameter that groups all constant terms of the first exponential factor in the equation for simple nucleation rate
- $\tau$ - Incubation time constant in the equation for simple nucleation rate
- $\Delta t_h$ - Nucleation hold time.

Dimensions (ellipsoidal semiaxes) of precipitate seeds

- a - semiaxis in the x-direction
- b - semiaxis in the y-direction
- c - semiaxis in the x-direction


## Appendix I: Example functions for $f_{\alpha}$, $f_{\beta}$, $f_{Landau}$, $H(\eta)$

$$
\begin{align}
f_{\alpha}(c_{\alpha}) = A_{2} c_{\alpha}^2 + A_{1} c_{\alpha} + A_{0}
\end{align}
$$

$$
\begin{align}
f_{\beta}(c_{\beta}) = B_{2} c_{\beta}^2 + B_{1} c_{\beta} + B_{0}
\end{align}
$$

$$
\begin{align}
f_{Landau}(\eta) = \eta^2  - 2\eta^3 +  \eta^4
\end{align}
$$

$$
\begin{align}
H(\eta) = 3 \eta^2 - 2 \eta^3
\end{align}
$$
