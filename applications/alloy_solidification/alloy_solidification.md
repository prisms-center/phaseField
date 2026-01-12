# PRISMS-PF: Alloy Solidification

This example application [1] implements a simple model to simulate the directional solidification of a binary alloy A-B in the dilute limit with component B acting as a solute in a matrix of A. 
The implemented model was introduced by Echebarria et al. [2] in 2004. In this model, latent heat is assumed to diffuse much faster than impurities and, therefore, the temperature field is considered
to be fixed by external conditions. In the default settings of the application, the simulation starts with a circular solid in the corner of an elongated system which evolves as a system is cooled under
a uniform thermal gradient and constant cooling rate. As this seed grows, three variables are tracked: an order parameter, $\phi$, that denotes whether the material is a liquid ($\phi=-1$) or solid ($\phi=1$),
a dimensionless supersaturation, $U$, and an auxiliary term, $\xi$. In addition, a solute concentration variable, $c$, is calculated during post-processing and outputted along with the aforementioned variables.

## Model
Consider a free energy expression given by:

$$
\begin{equation}
\mathcal{F}[\phi,c,T]= \int_{\Omega}  \left[\frac{\sigma}{2}|\nabla \phi|^2 +f(\phi, T_m) +  f_{AB}(\phi, c, T) \right]  ~dV,
\end{equation}
$$

where $\sigma$ is penalty coefficient for phase gradients. The term $f(\phi, T_m)$ is a symmetric double-well potential evaluated at the melting temperature of the pure material (A), $T_m$, and  the term $f_{AB}(\phi, c, T)$ accounts for the relative stability of the liquid and solid phases at different temperatures, according to the phase diagram. The double-well potential has the standard form with height $H$ given by

$$
\begin{equation}
f(\phi, T_m) = H(-\phi^2/2 + \phi^4/4).
\end{equation}
$$

For a dilute binary alloy, $f_{AB}$ can be written as

$$
\begin{equation}
f_{AB}(\phi, c, T) = f^A(T_m) - (T-T_m) s(\phi) +  \frac{R T_m}{v_o}(c\ln c -c) +\epsilon (\phi) c,
\end{equation}
$$

where $f^A(T_m)$ is the free energy density of pure A at its melting point, $v_0$ is the molar volume of A and $R$ is the gas constant. The functions $s(\phi)$ and $\epsilon(\phi)$ are interpolation functions for the entropy and internal energy of the solid and liquid phases, respectively. The general form of the coupled governing equations for the $\phi$ and $c$ is

$$
\begin{equation}
\frac{\partial \phi}{\partial t} = -K_\phi \frac{\delta \mathcal{F}}{\delta \phi}
\end{equation}
$$

and

$$
\begin{equation}
\frac{\partial c}{\partial t} = \nabla \cdot \left(\  M(\phi,c) \frac{\delta \mathcal{F}}{\delta c}  - \vec{\jmath}_{at} \right), 
\end{equation}
$$

where $K_\phi$ is a kinetic constant, $M(\phi,c)$ is the mobility of solute atoms and $\vec{\jmath}_{at}$ is a nonvariational anti-trapping solute current required to correct for spurious effects that arise from considering an interface thickness much larger than the physical solid-liquid interface. In the dilute limit, the solidus and liquidus lines of the $T$ vs. $c$ phase diagram are defined by the equations

$$
\begin{equation}
T_l=T_m-|m|c_l
\end{equation}
$$

and

$$
\begin{equation}
T_s=T_m-\frac{|m|}{k}c_s,
\end{equation}
$$

where $m$ is the liquidus slope and $k=c_s/c_l$ is the partition coefficient, which relates the equilibrium concentrations, $c_l$ and $c_s$, of the liquid and solid, respectively.

For the governing equations to simulate directional solidification, we follow the same approach of Ref. [2] by introducing a dimensionless supersaturation, $U$, instead of $c$. This supersaturation term is defined as

$$
\begin{equation}
U = \frac{e^u -1}{1-k},
\end{equation}
$$

where

$$
\begin{equation}
u(c,\phi) = \ln \left( \frac{2c}{c_l^0[1+k-(1-k)\phi]}\right).
\end{equation}
$$

The constant $c_l^0$, is the equilibrium liquidus concentration at reference temperature, $T_0$. We define it as $c_l^0  = c_\infty /k$, where $c_\infty$ is the concentration of the liquid far from the solid-liquid interface (equal to the average concentration of the alloy). Thus, the reference temperature is given by  $T_0=T_m - |m|c_l^0 =T_m-|m|c_\infty/k$.

## Governing Equations
After nondimensionalization (see Ref. [2] for derivation), the governing equations (in 2D) for $\phi$ and $U$ are given by

$$
\begin{equation}
\tau_\phi\frac{\partial  \phi}{\partial  t} = \xi(\phi,U)
\end{equation}
$$

and 

$$
\begin{equation}
\tau_U\frac{\partial  U}{\partial  t} = \nabla \cdot \left[ \tilde{D} \frac{1-\phi}{2} \nabla U - \vec{\jmath}_{at}^{\,U} \right] + \frac{1}{2}[1+(1-k)U]\frac{\partial \phi}{\partial t},
\end{equation}
$$

where

$$
\begin{align}
\tau_\phi=[1+(1-k) U ]a_s^2(\hat{n}),
\end{align}
$$

$$
\begin{align}
\tau_U=\frac{1+k}{2} - \frac{1-k}{2}\phi,
\end{align}
$$

$$
\begin{align}
\xi = & \nabla \cdot  \left( a_s^2(\hat{n}) \nabla \phi \right) +  \frac{\partial}{\partial x} \left[ |\nabla \phi|^2 a_s(\hat{n}) \frac{\partial a_s(\hat{n})}{\partial \left( \frac{\partial \phi}{\partial x} \right)} \right] +  \frac{\partial}{\partial y} \left[ |\nabla \phi|^2 a_s(\hat{n}) \frac{\partial a_s(\hat{n})}{\partial \left( \frac{\partial \phi}{\partial y} \right)} \right]\\
& +\phi-\phi^3 - \lambda(1-\phi^2)^2 \left[ U + U_\text{off} + \frac{\tilde{y} - \tilde{y}_0 - \tilde{V}_p t}{\tilde{l}_T} \right],
\end{align}
$$

and

$$
\begin{equation}
\vec{\jmath}_{at}^{\,U}=\frac{1}{2\sqrt{2}}[1+(1-k)U]\hat{n}\frac{\partial \phi}{\partial t}.
\end{equation}
$$

The function $a_s$ is the anisotropy factor for the solid-liquid interfacial energy, which depends on the outward normal (with respect to the solid) of the interface, 
$\hat{n}=-\nabla \phi / |\nabla \phi|$ (Note: The code uses the opposite sign convention). For a solid phase with $m$-fold symmetry this factor is given by 

$$
\begin{equation}
a_s(\hat{n})=1+\epsilon_m \cos[m(\theta-\theta_0)],
\end{equation}
$$

(In the implementation of the current model, $m$ is set to 4 and  $\theta_0=0$. For the purpose of computational efficiency, explicit calculation of trigonometric functions (and their inverse) is avoided. Thus, all sine and cosine terms with argument $m\theta$ are evaluated as $\sin(4\theta)=4\cos^3\theta\sin\theta-4\cos\theta\sin^3\theta$ and $\cos(4\theta)=\cos^4\theta -6\cos^2\theta\sin^2\theta-\sin^4\theta$, where $\sin\theta=\partial_y\phi / |\nabla \phi|$ and $\cos\theta=\partial_x\phi / |\nabla \phi|$.)

where $\epsilon_m$ determines the strength of the anisotropy, $\theta$ is the in-plane azimuthal angle of the normal vector with respect to the positive $x$-direction and $\theta_0$ is the reference orientation of the solid grains. The angle $\theta$ is related to the normal derivatives of $\phi$ at the interface via 

$$
\begin{equation}
\tan(\theta) = \frac{\partial \phi / \partial y}{\partial \phi / \partial x}. 
\end{equation}
$$

In 

$$
\begin{align}
\xi = & \nabla \cdot  \left( a_s^2(\hat{n}) \nabla \phi \right) +  \frac{\partial}{\partial x} \left[ |\nabla \phi|^2 a_s(\hat{n}) \frac{\partial a_s(\hat{n})}{\partial \left( \frac{\partial \phi}{\partial x} \right)} \right] +  \frac{\partial}{\partial y} \left[ |\nabla \phi|^2 a_s(\hat{n}) \frac{\partial a_s(\hat{n})}{\partial \left( \frac{\partial \phi}{\partial y} \right)} \right]\\
& +\phi-\phi^3 - \lambda(1-\phi^2)^2 \left[ U + U_\text{off} + \frac{\tilde{y} - \tilde{y}_0 - \tilde{V}_p t}{\tilde{l}_T} \right],
\end{align}
$$

$\lambda$ is a coupling constant defined as $\lambda=5\sqrt{2}W/(8d_0)$, where $W=(\sigma/H)^{1/2}$ is the equilibrium interface width and $d_0$ is the chemical capillary length, given by

$$
\begin{equation}
d_0=\frac{\gamma T_m}{L |m| (1-k) c_l^0}.
\end{equation}
$$

In the previous equation, $\gamma$ is the equilibrium surface tension and $L$ is the latent heat of fusion per volume. 
Finally, $\tilde{y}$, $\tilde{V}_p$, $\tilde{l}_T$ and $\tilde{D}$ are all dimensionless parameters, calculated by taking the unit length as $W$ and the unit time as $\tau_0=0.6267\lambda W^2/D$, 
where $D$ is the solute diffusivity in the liquid. The coordinate, $\tilde{y}$, represents the position along the direction of the thermal gradient, $\tilde{V}_p$ is the steady-state solidification speed, 
$\tilde{l}_T$ is the thermal length, calculated as

$$
\begin{align}
\tilde{l}_T=|m|(1-k)c_l^0/\tilde{G}
\end{align}
$$


where $\tilde{G}$ is the dimensionless thermal gradient, and $\tilde{D}$ is the dimensionless solute diffusivity in the liquid. Note that Eqs. 10 through 15 
are equivalent to Eqs. (132) and (133) from Ref. [2], except for the expression for the phase-field relaxation time $\tau_\phi$ which, for this application, was chosen to be $U$-dependent, 
as defined by Eq. (123) from Ref. [2].


Equation

$$
\begin{align}
\xi = & \nabla \cdot  \left( a_s^2(\hat{n}) \nabla \phi \right) +  \frac{\partial}{\partial x} \left[ |\nabla \phi|^2 a_s(\hat{n}) \frac{\partial a_s(\hat{n})}{\partial \left( \frac{\partial \phi}{\partial x} \right)} \right] +  \frac{\partial}{\partial y} \left[ |\nabla \phi|^2 a_s(\hat{n}) \frac{\partial a_s(\hat{n})}{\partial \left( \frac{\partial \phi}{\partial y} \right)} \right]\\
& +\phi-\phi^3 - \lambda(1-\phi^2)^2 \left[ U + U_\text{off} + \frac{\tilde{y} - \tilde{y}_0 - \tilde{V}_p t}{\tilde{l}_T} \right],
\end{align}
$$

can be simplified by explicitly writing $a_s(\hat{n})$ in terms of $\theta$.  We can evaluate the terms $\partial a_s(\theta)/\partial \left( \frac{\partial \phi}{\partial x} \right)$ and $\partial a_s(\theta)/\partial \left( \frac{\partial \phi}{\partial y} \right)$ by using the chain rule, i.e.,

$$
\begin{align}
\frac{\partial a_s(\theta)}{\partial \left( \frac{\partial \phi}{\partial x} \right)}=\frac{\partial a_s(\theta)}{\partial \theta} \frac{\partial \theta}{\partial \left( \frac{\partial \phi}{\partial x} \right)}\ \mathrm{and}\ \frac{\partial a_s(\theta)}{\partial \left( \frac{\partial \phi}{\partial y} \right)}=\frac{\partial a_s(\theta)}{\partial \theta} \frac{\partial \theta}{\partial \left( \frac{\partial \phi}{\partial y} \right)}
\end{align}
$$

along with 

$$
\begin{equation}
\tan(\theta) = \frac{\partial \phi / \partial y}{\partial \phi / \partial x}. 
\end{equation}
$$
See Appendix for details.

Also, the second and third terms on the right-hand side can be expressed using a divergence operator, allowing them to be grouped with the first term, which will simplify matters later. Carrying out these transformations yields:

$$
\begin{align}
\xi = & \nabla \cdot  \left[ \left(a_s^2(\theta) \frac{\partial \phi}{\partial x} + \epsilon_m m a_s(\theta) \sin \left[ m \left(\theta - \theta_0 \right) \right] \frac{\partial \phi}{\partial y}\right)\hat{x} \right.\\
& \left . + \left(a_s^2(\theta) \frac{\partial \phi}{\partial y} - \epsilon_m m a_s(\theta) \sin \left[ m \left(\theta - \theta_0 \right) \right] \frac{\partial \phi}{\partial x}\right)\hat{y}\right]\\
&+ \phi-\phi^3 - \lambda{(1-\phi^2)}^2 \left[ U + U_\text{off} + \frac{\tilde{y} - \tilde{y}_0 - \tilde{V}_p t}{\tilde{l}_T} \right].
\end{align}
$$

## Model Constants

- $\epsilon$: Strength of the anisotropy ($\epsilon_4$ for a solid with fourfold anisotropy)
- $k$: Partition coefficient
- $c_0$: Initial liquid concentration ($c_\infty$)
- $\lambda$: Coupling constant (setting this value fixes the interface width)
- $\tilde{D}$: Dimensionless solute diffusivity in the liquid phase.
- $\tilde{V}_p$: Dimensionless steady-state velocity of the tip.
- $\tilde{l}_T$: Dimensionless thermal length.
- $U_0$: Initial constitutional undercooling of the system ($U_0=-1$ sets the concentration of the liquid as $c_\infty$ and of the solid as $kc_\infty$)
- $U_{\text{off}}$: Undercooling offset that determines the initial temperature at the interface ($U_{\text{off}}=0$ sets it to the solidus temperature,  $U_{\text{off}}=1$ sets it to the liquidus temperature).
- $\tilde{y}_0$: Initial solid-liquid interface position relative to the bottom of the system ($\tilde{y}=0$)


## Time Discretization
Considering forward Euler explicit time stepping, we have the time-discretized kinetics equations:

$$
\begin{equation}
\phi^{n+1}=\phi^{n} + \frac{\xi^n}{\tau_\phi}\Delta t,
\end{equation}
$$

$$
\begin{align}
U^{n+1}=U^{n}+\frac{\Delta t}{\tau_U}\left[\nabla \cdot \left( \tilde{D}\frac{1-\phi^n}{2} \nabla U^n - \vec{\jmath}_{at}^{\ U} \right) + \frac{1}{2}[1+(1-k)U^n]\frac{\xi^n}{\tau_\phi} \right],
\end{align}
$$

and

$$
\begin{align}
\xi^{n+1} = & \nabla \cdot  \left[ \left(a_s^2(\theta^n) \frac{\partial \phi^n}{\partial x} + \epsilon_m m a_s(\theta^n) \sin \left[ m \left(\theta^n - \theta_0 \right) \right] \frac{\partial \phi^n}{\partial y}\right)\hat{x} \right.\\
& \left . + \left(a_s^2(\theta^n) \frac{\partial \phi^n}{\partial y} - \epsilon_m m a_s(\theta^n) \sin \left[ m \left(\theta^n - \theta_0 \right) \right] \frac{\partial \phi^n}{\partial x}\right)\hat{y}\right]\\
& +\phi^n-{(\phi^n)}^3 - \lambda {\left[1-{(\phi^n)}^2\right]}^2 \left[ U^n + U_\text{off} + \frac{\tilde{y} - \tilde{y}_0 - \tilde{V}_p t}{\tilde{l}_T} \right].
\end{align}
$$

## Weak Formulation
The weak formulation is obtained by multiplying the time-discretized equations by test function, $\omega$, and integrating over the volume, $\Omega$. For $\phi$ we get

$$
\begin{align}
\int_{\Omega}   \omega  \phi^{n+1}  ~dV = \int_{\Omega}   \omega \left(\phi^n + \frac{ \xi^n}{\tau_\phi}\Delta t\right) ~dV.
\end{align}
$$

$$
\begin{align}
r_{\phi} &= \left(\phi^n + \frac{ \xi^n}{\tau_\phi}\Delta t\right)
\end{align}
$$


For the weak form of 

$$
\begin{align}
U^{n+1}=U^{n}+\frac{\Delta t}{\tau_U}\left[\nabla \cdot \left( \tilde{D}\frac{1-\phi^n}{2} \nabla U^n - \vec{\jmath_{at}}^{\ U} \right) + \frac{1}{2}[1+(1-k)U^n]\frac{\xi^n}{\tau_\phi} \right],
\end{align}
$$

we employ the relation $\nabla \frac{1}{\tau_U}=\frac{1}{\tau_U^2}\frac{1-k}{2}\nabla\phi$ that results from substituting $\tau$ as defined by

$$
\begin{align}
\tau_U=\frac{1+k}{2} - \frac{1-k}{2}\phi,
\end{align}
$$

into the gradient of $1/\tau_U$:

$$
\begin{align}
\int_{\Omega}   \omega  U^{n+1}  ~dV =& 
\int_{\Omega} \omega \left( U^{n} + \frac{\Delta t}{2\tau_U\tau_\phi}[1+(1-k)U^n]\xi^n  - \frac{\Delta t (1-k)}{2\tau_U^2} \nabla \phi \cdot \left[\tilde{D}\frac{1-\phi^n}{2}\nabla U^n-\vec{\jmath}_{at}^{\,U}\right] \right) ~dV\\
&+\int_{\Omega}  \nabla  \omega  \cdot \left( -\frac{\Delta t}{\tau_U}\left[\tilde{D}(1-\phi^n)\nabla U^n-\vec{\jmath}_{at}^{\,U}\right] \right) ~dV.
\end{align}
$$

$$
\begin{align}
r_U &=  \left( U^{n} + \frac{\Delta t}{2\tau_U\tau_\phi}[1+(1-k)U^n]\xi^n  - \frac{\Delta t (1-k)}{2\tau_U^2} \nabla \phi \cdot \left[\tilde{D}\frac{1-\phi^n}{2}\nabla U^n-\vec{\jmath}_{at}^{\,U}\right] \right) 
\end{align}
$$

$$
\begin{align}
r_{Ux} &= \left( -\frac{\Delta t}{\tau_U}\left[\tilde{D}\frac{1-\phi^n}{2}\nabla U^n-\vec{\jmath}_{at}^{\,U}\right] \right)
\end{align}
$$

Finally, for $\xi$, we obtain

$$
\begin{equation}
\int_{\Omega}   \omega \xi^{n+1} ~dV =\int_{\Omega} \omega r_\xi ~dV + \int_{\Omega} \nabla \omega r_{\xi x} ~dV,
\end{equation}
$$

where

$$
\begin{equation}
r_\xi= \phi^n-(\phi^n)^3 - \lambda \left[1-(\phi^n)^2\right]^2 \left[ U^n + U_\text{off} + \frac{\tilde{y} - \tilde{y}_0 - \tilde{V}_p t}{\tilde{l}_T} \right]
\end{equation}
$$

$$
\begin{equation}
\begin{split}
r_{\xi x}= &-\left[a_s^2(\theta^n) \frac{\partial \phi^n}{\partial x} + \epsilon_m m a_s(\theta^n) \sin \left[ m \left(\theta^n - \theta_0 \right) \right] \frac{\partial \phi^n}{\partial y}\right]\hat{x}\\
&-\left[a_s^2(\theta^n) \frac{\partial \phi^n}{\partial y} - \epsilon_m m a_s(\theta^n) \sin \left[ m \left(\theta^n - \theta_0 \right) \right] \frac{\partial \phi^n}{\partial x}\right]\hat{y}
\end{split}
\end{equation}
$$


The above values of $r_{\phi}$, $r_{U}$, $r_{Ux}$,  $r_{\xi}$ and  $r_{\xi x}$ are used to define the residuals in the following parameters file:
`applications/alloy_solification/equations.cc`

## References
[1] Developed by Zhenjie Yao, Department of Material Science and Engineering, University of Michigan (2021).

[2] B. Echebarria, R. Folch, A. Karma, and M. Plapp, Quantitative phase-field model of alloy solidification, *Phys. Rev. E* **70**, 061604 (2004).


## Appendix

First, we express $\theta$ as a function of the gradient components:
$$
\begin{equation}
\theta = \arctan\left( \frac{\phi_y}{\phi_x} \right)
\end{equation}
$$

Applying the chain rule for $\theta$:
$$
\begin{align}
\frac{\partial \theta}{\partial \phi_x} &= \frac{1}{1 + (\phi_y/\phi_x)^2} \cdot \left( -\frac{\phi_y}{\phi_x^2} \right) = -\frac{\phi_y}{\phi_x^2 + \phi_y^2} = -\frac{\sin \theta}{|\nabla \phi|} \\
\frac{\partial \theta}{\partial \phi_y} &= \frac{1}{1 + (\phi_y/\phi_x)^2} \cdot \left( \frac{1}{\phi_x} \right) = \frac{\phi_x}{\phi_x^2 + \phi_y^2} = \frac{\cos \theta}{|\nabla \phi|}
\end{align}
$$

Now, applying the chain rule for $a(\theta)$:
$$
\begin{equation}
\frac{\partial a}{\partial \phi_x} = \frac{da}{d\theta} \frac{\partial \theta}{\partial \phi_x} = -\frac{a'(\theta) \sin \theta}{|\nabla \phi|}
\end{equation}
$$
$$
\begin{equation}
\frac{\partial a}{\partial \phi_y} = \frac{da}{d\theta} \frac{\partial \theta}{\partial \phi_y} = \frac{a'(\theta) \cos \theta}{|\nabla \phi|}
\end{equation}
$$
