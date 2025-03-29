# PRISMS-PF Application Formulation: grainGrowth

This example application implements a simple set of governing equations for isotropic grain growth. The model is a simplified version of the one in the following publication:
 Simulating recrystallization in titanium using the phase field method, S.P. Gentry and K. Thornton, *IOP Conf. Series: Materials Science and Engineering* 89 (2015) 012024.


Consider a free energy expression of the form:

$$
\begin{equation}
  \Pi(\eta_i, \nabla  \eta_i) = \int_{\Omega} \left[\sum_{i=1}^N \left(-\frac{1}{2}\eta_i^2+ \frac{1}{4}\eta_i^4 \right) + \alpha \sum_{i=1}^N \sum_{j>i}^N \eta_i^2 \eta_j^2 +\frac{1}{4} \right] +  \frac{\kappa}{2} \sum_{i=1}^N |\nabla \eta_i|^2    ~dV 
\end{equation}
$$

where $\eta_i$ is one of $N$ structural order parameters, $\alpha$ is the grain interaction coefficient, and $\kappa$ is the gradient energy coefficient.
	
## Variational treatment
The driving force for grain evolution is determined by the variational derivative of the total energy with respect to each order parameter:

$$
\begin{equation}
\mu = \frac{\delta \Pi}{\delta \eta_i} = \left( -\eta_i + \eta_i^3 + 2 \alpha \eta_i \sum_{j \ne i}^N \eta_j^2 - \kappa \nabla^2 \eta_i \right)
\end{equation}
$$

## Kinetics
The order parameter for each grain is unconserved, and thus their evolution can be described by Allen-Cahn equations:

$$
\begin{equation}
\frac{\partial \eta_i}{\partial t} = -L \mu = \left( -\eta_i + \eta_i^3 + 2 \alpha \eta_i \sum_{j \ne i}^N \eta_j^2 - \kappa \nabla^2 \eta_i \right)
\end{equation}
$$

where $L$ is the constant mobility. 

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
