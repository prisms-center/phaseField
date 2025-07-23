# PRISMS-PF Application Formulation: srx_with_nucleation

This example expands on the simple isotropic grain growth example (applications/grainGrowth) by introducing a stored energy driving force and nucleation of unstrained grains.
This model is based on that of Huang et al., "Phase-field modeling of stored-energy-driven grain growth with intra-granular variation in dislocation density," *Modelling Simul. Mater. Sci. Eng.* **32** 045011 (2024), which is in turn based on that of Susan P. Gentry and Katsuyo Thornton, "Sensitivity analysis of a phase field model for static recrystallization of deformed microstructures," *Modelling Simul. Mater. Sci. Eng.* **28** 065002 (2020).


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

Each grain is assigned a dislocation density $\rho_i$. By default, the per-grain dislocation density fields do not change with time; only the overall interpolated dislocation density changes as the order parameters evolve. Static recovery can be simulated by making the dislocation density decrease with time, if a recovery rate law is defined. To do this, an evolution equation for the dislocation density fields should be defined in `equations.cc`.
	
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

## Time discretization
Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:

$$
\begin{equation}
 \eta_i^{n+1} = \eta_i^{n} - \Delta t \frac{\partial \eta_i}{\partial t}
\end{equation}
$$
 
## Weak formulation
In the weak formulation, considering an arbitrary variation $w$, the above equation can be expressed as a residual equation:

$$
\begin{align}
\int_{\Omega} w \eta_i^{n+1} dV &= \int_{\Omega}\left[ w \eta_i^{n} - w \Delta t L \left(m_0\left(-\eta_i + \eta_i^3 + 2 \alpha \eta_i \sum_{j \ne i}^N \eta_j^2 \right) - \kappa \nabla^2 \eta_i
                                     + Gb^2 \frac{\eta_i}{\sum_{j=1}^{N}} \left(\rho_i - \rho\right) \right) \right]dV \\
&= \int_{\Omega} \left[ w \left( \eta_i^{n} - \Delta t L m_0\left(-\eta_i + \eta_i^3 + 2 \alpha \eta_i \sum_{j \ne i}^N \eta_j^2 \right) + \Delta t L Gb^2 \frac{\eta_i}{\sum_{j=1}^{N}} \left(\rho_i - \rho\right) \right) + \nabla w (-\Delta t L \kappa) \cdot (\nabla \eta_i^{n}) \right]dV \quad [\kappa \nabla \eta_i \cdot n = 0 \text{on} \partial \Omega]
\end{align}
$$

$$
\begin{align}
r_{\eta_i} &= \eta_i^{n} - \Delta t L m_0\left(-\eta_i + \eta_i^3 + 2 \alpha \eta_i \sum_{j \ne i}^N \eta_j^2 \right) + \Delta t L Gb^2 \frac{\eta_i}{\sum_{j=1}^{N}} \left(\rho_i - \rho\right) \\
r_{\eta_i x} &= (-\Delta t L \kappa) \cdot (\nabla \eta_i^{n})
\end{align}
$$

The above values of  $r_{\eta_i}$ and $r_{\eta_i x}$ are used to define the residuals in the following C++ file: 
`applications/srx_with_nucleation/equations.cc`
