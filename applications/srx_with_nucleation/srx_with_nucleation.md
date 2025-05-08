# PRISMS-PF Application Formulation: srx_with_nucleation

This example expands on the simple isotropic grain growth example (applications/grainGrowth) by introducing a stored energy driving force and nucleation of unstrained grains.
This model was used for simulation of static recrystallization in the publication:
 Susan P. Gentry and Katsuyo Thornton, "Sensitivity analysis of a phase field model for static recrystallization of deformed microstructures," *Modelling Simul. Mater. Sci. Eng.* **28** 065002 (2020).


Consider a free energy expression of the form:

$$
\begin{equation}
  F(\eta_i, \nabla  \eta_i) = \int_{\Omega} \left[m_0\left[\sum_{i=1}^{N} \left(-\frac{1}{2}\eta_i^2 + \frac{1}{4}\eta_i^4 \right)
                            + \alpha \sum_{i=1}^{N} \sum_{j>i}^{N}\left( \eta_i^2 \eta_j^2 \right) + \frac{1}{4} \right]
							+ \frac{\kappa}{2} \sum_{i=1}^N \left|\nabla \eta_i\right|^2  + f_s\right]  dV 
\end{equation}
$$

where $\eta_i$ is one of $N$ structural order parameters, $\alpha$ is the grain interaction coefficient, and $\kappa$ is the gradient energy coefficient.
$f_s$ is a contribution from stored energy due to the presence of dislocations in the deformed grain structure.

The contribution from stored energy, $f_s$, is:
$$
\begin{equation}
  f_s = \frac{1}{2}\rho G b^2
\end{equation}
$$
where $G$ is the shear modulus, $b$ is the Burger's vector, and $\rho$ is the dislocation density interpolated based on the order parameters:
$$
\begin{equation}
  \rho(\eta_i) = \frac{\sum_{i=1}^{N} \eta_i^2\rho_i}{\sum_{j=1}^{N} \eta_j^2}
\end{equation}
$$
Each grain is assigned a dislocation density $\rho_i$, which is assumed constant by default.
Static recovery can be simulated by making the dislocation density decrease with time, if a recovery rate law is defined. 
	
## Variational treatment
The driving force for grain evolution is determined by the variational derivative of the total energy with respect to each order parameter:

$$
\begin{equation}
\mu = \frac{\delta F}{\delta \eta_i} = m_0\left(-\eta_i + \eta_i^3 + 2 \alpha \eta_i \sum_{j \ne i}^N \eta_j^2 \right) - \kappa \nabla^2 \eta_i
                                     + Gb^2 \frac{\eta_i}{\sum_{j=1}^{N}} \left(\rho_i - \rho\right)
\end{equation}
$$

## Kinetics
The order parameters are unconserved, so their evolution is governed by a series of Allen-Cahn equations with constant mobility $L$:
$$
\begin{equation}
\frac{\partial \eta_i}{\partial t} = -L \mu = -L\left[m_0\left(-\eta_i + \eta_i^3 + 2 \alpha \eta_i \sum_{j \ne i}^N \eta_j^2 \right) - \kappa \nabla^2 \eta_i
                                     + Gb^2 \frac{\eta_i}{\sum_{j=1}^{N}} \left(\rho_i - \rho\right)\right]
\end{equation}
$$

#TODO finish the definition of the weak form

## Time discretization
Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:

$$
\begin{align}
 \eta_i^{n+1} &= \eta_i^{n} - \Delta t L~\left( -\eta_i^n + (\eta_i^n)^3 + 2 \alpha \eta_i^n \sum_{j \ne i}^N (\eta^n_j)^2 - \kappa \nabla^2 \eta^n_i \right)
\end{align}
$$
 
## Weak formulation
In the weak formulation, considering an arbitrary variation $w$, the above equation can be expressed as a residual equation:

$$
\begin{align}
\int_{\Omega}   w \eta_i^{n+1} ~dV&= \int_{\Omega}   w \eta_i^{n} - w \Delta t L~\left( -\eta_i^n + (\eta_i^n)^3 + 2 \alpha \eta_i^n \sum_{j \ne i}^N (\eta^n_j)^2 - \kappa \nabla^2 \eta^n_i \right) ~dV
\end{align}
$$

$$
\begin{align}
&= \int_{\Omega}   w ( \eta^{n} - \Delta t L~\left( -\eta_i^n + (\eta_i^n)^3 + 2 \alpha \eta_i^n \sum_{j \ne i}^N (\eta^n_j)^2\right) + \nabla w (-\Delta t L \kappa)~ \cdot (\nabla \eta_i^{n}) ~dV \quad [\kappa \nabla \eta_i \cdot n = 0 ~ \text{on} ~ \partial \Omega]
\end{align}
$$

$$
\begin{align}
r_{\eta_i} &= \eta^{n} - \Delta t L~\left( -\eta_i^n + (\eta_i^n)^3 + 2 \alpha \eta_i^n \sum_{j \ne i}^N (\eta^n_j)^2\right)
\end{align}
$$

$$
\begin{align}
r_{\eta_i x} &=  (-\Delta t L \kappa)~ \cdot (\nabla \eta_i^{n})
\end{align}
$$

The above values of  $r_{\eta_i}$ and $r_{\eta_i x}$ are used to define the residuals in the following parameters file: 
`applications/grainGrowth/equations.h`
