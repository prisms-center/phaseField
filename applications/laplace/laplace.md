# PRISMS-PF: Laplace's Equation

This example application solves [Laplace's equation](https://en.wikipedia.org/wiki/Laplace%27s_equation) with Dirichlet boundary values.

$$
\begin{align}
\nabla\cdot\nabla u_{\Omega} = 0 \\
u_{\partial \Omega} = \Gamma_D
\end{align}
$$

where $\Gamma_D$ are the Dirichlet boundary values.

## Governing equations
In the weak form (and multiplying by -1), we have

$$
\begin{align}
\int_{\Omega} \nabla w \cdot \nabla u ~dV = \int_{\Omega}0~dV
\end{align}
$$



## Expressions
In PRISMS-PF, linear equations are solved provided expressions for the left and right hand side of equations of the form $Ax=b$. From the weak form, we find

$$
\begin{align}
r_{ux}^{LHS} = \nabla u^*
\end{align}
$$

$$
\begin{align}
r_{u}^{RHS} = 0
\end{align}
$$

The ${}^*$ indicates a trial solution.

The above equations are used to define expressions in the file:
`applications/laplace/custom_pde.h`
