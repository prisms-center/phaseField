# PRISMS-PF Application Formulation: fickianDiffusion

In this example application, we implement a Fick's Law for a single component. Two time-dependent Gaussian source terms add concentration over the course of the simulation.

## Kinetics
The Parabolic PDE for diffusion is given by:

$$
\begin{align}
  \frac{\partial c}{\partial t} &= -\nabla \cdot (-D~\nabla c) + f
\end{align}
$$

where $D$ is the diffusion constant and $f$ is a source term. In this example, $f$ is given by a pair of Gaussian expressions:

$$
\begin{multline}
f = A_1 \exp\left(-\left(\frac{t-t_1}{\tau_1}\right)^2\right)\exp\left(-\left(\frac{x-x_1}{L_1}\right)^2-\left(\frac{y-y_1}{L_1}\right)^2\right) + A_2 \exp\left(-\left(\frac{t-t_2}{\tau_2}\right)^2\right)\exp\left(-\left(\frac{x-x_2}{L_2}\right)^2-\left(\frac{y-y_2}{L_2}\right)^2\right)
 \end{multline}
$$

 where  $A_1$, $_2$, $t_1$, $t_2$, $\tau_1$, $\tau_2$, $x_1$, $x_2$, $y_1$, $y_2$, $L_1$, and $L_2$ are constants.

## Time discretization
Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:

$$
\begin{align}
 c^{n+1} &= c^{n} + (\Delta t D)~\Delta c^n + \Delta t f^n
\end{align}
$$

## Weak formulation
In the weak formulation, considering an arbitrary variation $w$, the above equation can be expressed as a residual equation:

$$
\begin{align}
\int_{\Omega}   w c^{n+1} ~dV &= \int_{\Omega}   w c^{n} + w (\Delta t D) \Delta c^n + w \Delta t f^n ~dV
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}   w (c^{n} + \Delta t f^n) - \nabla w  \cdot (\Delta t D) \nabla c^n ~dV + \int_{\partial \Omega}   w  (\Delta t D) \nabla c^n \cdot n ~dS
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}   w (c^{n} + \Delta t f^n) - \nabla w  \cdot (\Delta t D) \nabla c^n ~dV + \int_{\partial \Omega}   w  (\Delta t D) j^n  ~dS
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}   w (c^{n}+\Delta t f^n) + \nabla w  \cdot (-\Delta t D) \nabla c^n ~dV \quad [\text {assuming flux}~j=0 ]
\end{align}
$$

$$
\begin{align}
r_c &= c^{n}+\Delta t f^n
\end{align}
$$

$$
\begin{align}
r_{cx} &= (-\Delta t D) \nabla c^n
\end{align}
$$

The above values of  $r_{c}$ and $r_{c x}$ are used to define the residuals in the following parameters file:
`applications/fickianFlux/parameters.h`
