# PRISMS-PF Application Formulation: alloySolidification_uniform

This example application implements a simple model to simulate solidification of a binary alloy A-B in the dilute limit with component B acting as a solute in a matrix of A. The implemented model was introduced by Karma [1] in 2001. In this model, latent heat is assumed to diffuse much faster than impurities and, therefore, the temperature field is considered to be fixed by external conditions. In contrast to *alloySolidification*, this application considers solidification under uniform temperature and no diffusion in the solid.  In the default settings of the application, the simulation starts with a circular solid in the corner of a square system. The evolution of the system is calculated for an initial scaled supersaturation value, $\Omega$. As this seed grows, three variables are tracked, an order parameter, $\phi$, that denotes whether the material a liquid ($\phi=-1$) or solid ($\phi=1$), the solute concentration, $c$, and an auxiliary term, $\xi$. 

## Model
The free energy of the system can be written as the functional [2]

$$
\begin{equation}
\mathcal{F}[\phi,c,T]= \int_{\Omega}  \left[\frac{\sigma}{2}|\nabla \phi|^2 +Hf(\phi, T_M) +  f_{AB}(\phi, c, T) \right]  ~dV,
\end{equation}
$$


where $\sigma$ is penalty coefficient for phase gradients. The term $f(\phi, T_M)$ is a symmetric double-well term evaluated at the melting temperature of the pure material $A$,  $T_M$. The constant $H$ is the height of the well and $f_{AB}(\phi, c, T)$ accounts for the relative stability of the liquid and solid phases at different temperatures, according to the phase diagram. The double-well term has the standard form given by

$$
\begin{equation}
f(\phi, T_M) = -\phi^2/2 + \phi^4/4.
\end{equation}
$$

For a dilute binary alloy, $f_{AB}$ can be written as

$$
\begin{equation}
f_{AB}(\phi, c, T) = f^A(T_M) - (T-T_M) s(\phi) +  \frac{R T_M}{v_o}(c\ln c -c) +\epsilon (\phi) c,
\end{equation}
$$

where $f^A(T_M)$ is the free energy density of pure A at its melting point, $v_0$ is the molar volume of A and $R$ is the gas constant. The functions $s(\phi)$ and $\epsilon (\phi)$ are interpolation functions for the entropy and internal energy of the solid and liquid phases, respectively. The general form of the coupled governing equations for the $\phi$ and $c$ is

$$
\begin{equation}
\frac{\partial \phi}{\partial t} = -K_\phi \frac{\delta \mathcal{F}}{\delta \phi}
\end{equation}
$$

and

$$
\begin{equation}
\frac{\partial c}{\partial t} = \nabla \cdot \left(\  M(\phi,c) \frac{\delta \mathcal{F}}{\delta c}  - \vec{j}_{at} \right), 
\end{equation}
$$

where $K_\phi$ is a kinetic constant, $M(\phi,c)$ is the mobility of solute atoms and $\vec{j}_{at}$ is a nonvariational anti-trapping solute current required to correct for spurious effects that arise from considering an interface thickness much larger than the physical solid-liquid interface. in the dilute limit, the solidus and liquidus lines of the $T$ vs $c$ phase diagram are defined by the equations

$$
\begin{equation}
T_l=T_M-|m|c_l
\end{equation}
$$

and

$$
\begin{equation}
T_s=T_M-\frac{|m|}{k}c_s,
\end{equation}
$$

where $m$ is the liquidus slope and $k=c_s/c_l$ is the partition coefficient, which relates the equilibrium concentrations, $c_l$ and $c_s$, of the liquid and solid, respectively.

Considering solidification at a temperature $T_0<T_M$, we define $c_l^0=c(T_0)$ as the equilibrium liquid concentration and $kc_l^0$ as the equilibrium solid concentration. After nondimensionalization, the coupled governing equations for the $\phi$ and $c$ can be written as [1]

$$
\begin{equation}
\tau\frac{\partial \phi}{\partial t} = \xi(\phi,c),
\end{equation}
$$

where

$$
\begin{equation}
\xi(\phi,c) = -\frac{\partial f}{\partial \phi} - \frac{\lambda}{1-k} g'(\phi)(e^u - 1) + W^2\nabla^2 \phi 
\end{equation}
$$

and

$$
\begin{equation}
\frac{\partial c}{\partial t} = \nabla \cdot \vec{j},
\end{equation}
$$

where

$$
\begin{equation}
\vec{j}=-D c q(\phi)\nabla u - aWc_l^0(1-k)e^u\frac{\partial \phi}{\partial t}\frac{\nabla \phi}{|\nabla \phi|},
\end{equation}
$$

and

$$
\begin{equation}
u(c,\phi) = \ln \left( \frac{2c/c_l^0}{1+k-(1-k)h(\phi)}\right).
\end{equation}
$$

The constant  $\lambda$ is defined as $\lambda=a_1W/d_0$, where $d_0$ is the microscopic capillary length. The constants $W$ and $\tau$ are the unit length and time, respectively and $D$ is the diffusivity in the liquid. The value of $D$ is set to $D=a_2W^2\lambda /\tau$ so that interface kinetics are eliminated from the velocity-dependent Gibbs-Thompson condition. The values of $a$, $a_1$, $a_2$ are given in the Model Constants section. The functions $g'(\phi)$, $h(\phi)$ and $q(\phi)$ are given by

$$
\begin{equation}
g'(\phi)=(1-\phi^2)^2, 
\end{equation}
$$

$$
\begin{equation}
q(\phi)=\frac{1-\phi}{1+k-(1-k)h(\phi)},
\end{equation}
$$

and

$$
\begin{equation}
h(\phi) = \phi.
\end{equation}
$$

Crystalline anisotropy is introduced by generalizing Eqs.

$$
\begin{equation}
\tau\frac{\partial \phi}{\partial t} = \xi(\phi,c),
\end{equation}
$$

and 

$$
\begin{equation}
\xi(\phi,c) = -\frac{\partial f}{\partial \phi} - \frac{\lambda}{1-k} g'(\phi)(e^u - 1) + W^2\nabla^2 \phi 
\end{equation}
$$

to

$$
\begin{equation}
\tau(\theta)\frac{\partial \phi}{\partial t} = \xi(\phi,c),
\end{equation}
$$

where

$$
\begin{equation}
\xi(\phi,c) = -f'(\phi) - \frac{\lambda}{1-k} g'(\phi)(e^u - 1) + \nabla \cdot [W(\theta)^2\nabla \phi]-\frac{\partial}{\partial x} \left[W(\theta)W'(\theta)\frac{\partial \phi}{\partial y}\right] + \frac{\partial}{\partial y} \left[W(\theta)W'(\theta)\frac{\partial \phi}{\partial x}\right],
\end{equation}
$$

where $\theta$ is the angle between the outward normal of the solid-liquid interface and the positive $x$ axis. The constants $W$ and $\tau$ now depend on this angle and are given by $W(\theta)=Wa_s(\theta)$ and $\tau(\theta)=\tau a_s(\theta)^2$, respectively, where $a_s(\theta)=1+\epsilon_4 \cos(4\theta)$.

The initial condition is set by placing a solid seed in an undercooled system with a uniform scaled supersaturation value, $\Omega=(c_l^0-c_\infty)/[c_l^0(1-k)]$.

## Model Constants

| Symbol | Value | Description |
|--------|-------|-------------|
|$c_l^0$ | 1     | Reference concentration (equilibrium liquid concentration at $T_0$) |
|$k$     | 0.15  | Partition coefficient |
| $W$    | 1     | Unit length |
| $\tau$ | 1     | Unit time |
| $d_0/W$| 0.277 |Microscopic capillary length (with respect to W) |
| $a$ | $1/(2\sqrt{2})$ |Antitraping term constant|
|$a_1$| 0.8839 | Coefficient $a_1$ | 
| $a_2$ | 0.6267 | Coefficient $a_2$ |
| $\lambda$ | $a_1W/d_0$ | Parameter $\lambda$ |
| $D$ | $a_2\lambda W^2/\tau$ | Solute diffusivity in the liquid |
| $\epsilon_4$ | 0.02 | Anisotropy strength|
| $\Omega$ | 0.55 | Scaled supersaturation|


## Time Discretization
Considering forward Euler explicit time stepping, we have the time discretized kinetics equations:

$$
\begin{equation}
\phi^{n+1}=\phi^{n} + \frac{\xi^n}{\tau^n}\Delta t,
\end{equation}
$$

$$
\begin{equation}
c^{n+1}=c^{n}+\Delta t \left \( D c^n q(\phi^n)\nabla(u^n) + aWc_l^0(1-k)(e^u)^n \left(\frac{\partial \phi}{\partial t}\right)^n\frac{\nabla \phi^n}{|\nabla \phi^n|} \right \) 
\end{equation}
$$

and

$$
\begin{equation}
\begin{split}
\xi(\phi,c) = & -f'(\phi^n) - \frac{\lambda}{1-k} g'(\phi^n)[(e^u)^n - 1]\\
& + \nabla \cdot [W(\theta^n)^2\nabla \phi^n]-\frac{\partial}{\partial x} \left[W(\theta^n)W'(\theta^n)\frac{\partial \phi^n}{\partial y}\right] + \frac{\partial}{\partial y} \left[W(\theta^n)W'(\theta^n)\frac{\partial \phi^n}{\partial x}\right].
\end{split}
\end{equation}
$$

## Weak Formulation
The weak form of the time-discretized equations for $\phi$,  $c$, and, $\xi$ is

$$
\begin{equation}
\int_{\Omega}   \omega  \phi^{n+1}  ~dV = \int_{\Omega}   \omega \left(\phi^n + \frac{ \xi^n}{\tau(\theta^n)}\Delta t\right) ~dV,
\end{equation}
$$

$$
\begin{align}
\int_{\Omega}   \omega  c^{n+1}  ~dV =& 
\int_{\Omega} \omega  c^{n} ~dV
\end{align}
$$

$$
\begin{align}
r_{\phi} &= \left(\phi^n + \frac{ \xi^n}{\tau(\theta^n)}\Delta t\right)
\end{align}
$$

$$
\begin{align}
r_c &= c^{n} 
\end{align}
$$

$$
\begin{align}
&+\int_{\Omega}  \nabla  \omega  \cdot \left[-\Delta t D\left(q(\phi^n)\nabla c^n + \frac{(1-k)q(\phi^n)c^n\nabla(\phi^n)}{1+k-(1-k)\phi^n}\right)-\Delta t aWc_l^0(1-k)(e^u)^n \left(\frac{\partial \phi}{\partial t}\right)^n\frac{\nabla \phi^n}{|\nabla \phi^n|} \right] ~dV,
\end{align}
$$

$$
\begin{align}
r_cx &= \left[-\Delta t D\left(q(\phi^n)\nabla c^n + \frac{(1-k)q(\phi^n)c^n\nabla(\phi^n)}{1+k-(1-k)\phi^n}\right)-\Delta t aWc_l^0(1-k)(e^u)^n \left(\frac{\partial \phi}{\partial t}\right)^n\frac{\nabla \phi^n}{|\nabla \phi^n|} \right]
\end{align}
$$

and

$$
\begin{align}
\int_{\Omega}   \omega \xi^{n+1} ~dV =\int_{\Omega} \omega r_\xi ~dV + \int_{\Omega} \nabla \omega r_{\xi x} ~dV
\end{align}
$$

where

$$
\begin{equation}
r_\xi= -f'(\phi^n) - \frac{\lambda}{1-k} g'(\phi^n)[(e^u)^n - 1]
\end{equation}
$$

and 

$$
\begin{align}
r_{\xi x}= &-\left[W(\theta^n)^2\frac{\partial \phi^n}{\partial x}-W(\theta^n)W'(\theta^n)\frac{\partial \phi^n}{\partial y}\right]\hat{x}
\end{align}
$$

$$
\begin{align}
&-\left[W(\theta^n)^2\frac{\partial \phi^n}{\partial y}+W(\theta^n)W'(\theta^n)\frac{\partial \phi^n}{\partial x}\right]\hat{y}
\end{align}
$$


The above values of $r_{\phi}$, $r_{c}$, $r_{cx}$,  $r_{\xi}$ and  $r_{\xi x}$ are used to define the residuals in the following parameters file:
`applications/alloySolification_uniformT/equations.cc`

## References
[1] A. Karma, Phase-Field Formulation for Quantitative Modeling of Alloy Solidification, *Phys. Rev. Lett.*  **87**, 115701 (2001).

[2] B. Echabarria, R. Folch, A. Karma, and M. Plapp, Quantitative phase-field model of alloy solidification, *Phys. Rev. E* **70**, 061904 (2004).
