Analytical solution in lecture (lecture 7?)

Q=0

![[Reactor Physics/numerical exercise/code/ReactorSimulation.jl/docs/foreign/L07_Diffusion_part2 (1).pdf#page=15&rect=4,15,716,516|L07_Diffusion_part2 (1), p.15]]

- [ ] How to deal with half point at the boundary condition?
- [ ] Flux is 0 at the extrapolated length
- [ ] Length of slab? ![[Pasted image 20241008135044.png]]


![[Pasted image 20241008140047.png]]

The boundary condition is from Lecture 6, $J^+$ / $J^-$

# Documentation
ODE:
$$
D \frac{ d^2 }{ dx^2 } \Phi  -  \Sigma_{a}\Phi = 0
$$
Left boundary
$$
D \frac{ d^2 }{ dx^2 } \Phi(0)  -  \Sigma_{a}\Phi(0) = \frac{S}{2}
$$
Right Boundary
$$
J^{+}(0) = \frac{\phi(0)}{4} - D \left. \frac{d\phi(x)}{dx} \right|_{x=0} = 0
$$

[[Implementing the Boundary Conditions (outdated)]]

# Davids Solution
$$
\vec{S} = \begin{pmatrix}
- \frac{S}{dx \cdot 2D} \\
0 \\
\dots \\
0
\end{pmatrix}
$$
%% `boundary = 1/dx^2 * [1; zeros(n-2); (2*D/((1/4+D/dx)*dx)-1)] `%%
$$
\vec{B} = \frac{1}{dx^2} \begin{pmatrix}
1  \\
0 \\
\dots \\
\frac{2D \cdot dx}{\frac{1}{4}+\frac{D}{dx}}-1
\end{pmatrix}
$$


# Template Numerical Exercise 1
Relationship between $\mathrm{\Phi}_i,\ \mathrm{\Phi}_{i+1},\ \mathrm{\Phi}_{i-1}$ at any point within the material:
$$
\frac{1}{dx^2}(\Phi_{i-1} - 2 \Phi_{i} + \Phi_{i+1}) - \frac{1}{L} \Phi_{i} = 0
$$  

Coefficients of the matrix A (4 significant digits), for a mesh size of 0.1cm:\\

Coef ${} A_{i,i} = -200.2411 cm^{-2}$: \\

Coef $A_{i-1,i} = 99.9999 cm^{-2}$: \\

Coef $A_{i+1,i} = 99.9999 cm^{-2}$: \\

```
A[5, 5] = -200.24119999999996 cm^-2
A[5, 6] = 99.99999999999999 cm^-2
A[5, 4] = 99.99999999999999 cm^-2
A[4, 5] = 99.99999999999999 cm^-2
A[6, 5] = 99.99999999999999 cm^-2
```

# Correct Boundary Condition
Right hand side
$$
\frac{1}{dx^2} ( D \phi_{n-1} - (D + dx BC_{n} + \Sigma_{a} dx^2)\phi_{n}) = 0
$$
Left hand side:
$$
\frac{1}{dx^2}(-(D + \Sigma_{a} dx^2) \phi_{1} + D\phi_{2}) = -\frac{S}{2dx}
$$
Middle:
$$
\frac{D}{dx^2}(\Phi_{i-1} - 2 \Phi_{i} + \Phi_{i+1}) - \Sigma_{a} \Phi_{i} = 0
$$
