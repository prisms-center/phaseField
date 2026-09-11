For an isotropic material with $E = 1.0$ and $\nu = 0.25$, the **plane stress** stiffness matrix is
$$\mathbf{Q} =  \begin{bmatrix}  1.0667 & 0.2667 & 0 \\  0.2667 & 1.0667 & 0 \\  0 & 0 & 0.4000  \end{bmatrix}$$
and the **plane strain** stiffness matrix is
$$\mathbf{C} =  \begin{bmatrix}  1.2 & 0.4 & 0.4 & 0 \\  0.4 & 1.2 & 0.4 & 0 \\  0.4 & 0.4 & 1.2 & 0 \\  0 & 0 & 0 & 0.4  \end{bmatrix}$$

Applying a uniform uniaxial strain $\varepsilon_{11} = 0.01$ results in the following stress calculations:

- **Plane Stress:**
$$\begin{bmatrix}  \sigma_{11} \\ \sigma_{22} \\ \sigma_{12}  \end{bmatrix}  =  \begin{bmatrix}  1.0667 & 0.2667 & 0 \\  0.2667 & 1.0667 & 0 \\  0 & 0 & 0.4000  \end{bmatrix} \begin{bmatrix}  0.01 \\ 0 \\ 0  \end{bmatrix} = \begin{bmatrix}  0.010667 \\ 0.002667 \\ 0  \end{bmatrix}$$

- **Plane Strain:**
$$\begin{bmatrix}  \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{12}  \end{bmatrix}  =  \begin{bmatrix}  1.2 & 0.4 & 0.4 & 0 \\  0.4 & 1.2 & 0.4 & 0 \\  0.4 & 0.4 & 1.2 & 0 \\  0 & 0 & 0 & 0.4  \end{bmatrix} \begin{bmatrix}  0.01 \\ 0 \\ 0 \\ 0  \end{bmatrix} = \begin{bmatrix}  0.012 \\ 0.004 \\ 0.004 \\ 0  \end{bmatrix}$$

These analytical values match the output of the code.
