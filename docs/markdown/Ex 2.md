$$
D \frac{d^2}{dx^2}\phi  - \Sigma_{a}\phi = 0
$$

# Boundary Conditions
left
$$
J^{-}(0) = \frac{\phi(0)}{4} + D \left. \frac{d\phi(x)}{dx} \right|_{x=0} = 0
$$
right
$$
J^{+}(0) = \frac{\phi(a)}{4} - D \left. \frac{d\phi(x)}{dx} \right|_{x=a} = 0
$$

# Bare Reactor Analytical Solution
For the bare reactor the Equation becomes
$$
\frac{d^2 \phi}{d z^2}+B_m^2 \phi=0
$$
where $B_{m}$ is the material Buckling. 
It is solved by $\phi(z)=A \cos (B z)+C \sin (B z)$. Due to the symmetric boundary conditions of the problem, we know that $C=0$. The geometric buckling is $B^2=\left(\frac{\pi}{H+2 d}\right)^2$.
A steady state is reached, when the material Buckling is equal to the numerical buckling, in this case, $k_{eff}$ can be calculated by
$$
k_{\mathrm{eff}}=\frac{k_{\infty}}{1+L^2 B^2}
$$

# Bare Reactor Numerical Solution
```
JL = (P[2] - P[1]) / dx = 0.12542781405179504
JR = (P[end - 1] - P[end]) / dx = 0.12542781405179393
```

$$
J = \frac{\phi_{2}-\phi_{1}}{dx}= 0.12543
$$
# Reflected Reactor Numerical Solution



# Boundary Condition on the Sides
[[TP_I_intro_slides (new).pdf#page=10&rect=91,375,409,443|TP_I_intro_slides (new), p.10]]
Current from the right side is 0
$$
J_{x_{i+1 / 2}}^{-}=\frac{\Phi_{i+1 / 2}}{4}+D_i \frac{\Phi_{i+1 / 2}-\Phi_i}{\Delta x}=0
$$
reanranging to $\phi_{i+\frac{1}{2}}$:
$$
\Phi_{i+1 / 2}=\frac{1}{\left(\frac{1}{4}+\frac{D_i}{\Delta x}\right)} \frac{D_i}{\Delta x} \Phi_i=\frac{1}{\left(\frac{\Delta x}{4 D_i}+1\right)} \Phi_i
$$
plugging into the formula for the derivative:
$$
D_i \frac{\Phi_{i+1 / 2}-\Phi_i}{\Delta x / 2}=\frac{D_{i, j}}{\Delta x / 2} \Phi_i\left(\frac{1}{\left(\frac{\Delta x}{4 D_i}+1\right)}-1\right)=-\frac{1}{2\left(\frac{\Delta x}{4 D_i}+1\right)} \Phi_i
$$
let $B C_i=\frac{1}{2\left(\frac{\Delta x}{4 D_i}+1\right)}$.
and as before the backwards derivative is:
$$
-D_i \frac{\Phi_i-\Phi_{i-\frac{1}{2}}}{\Delta x / 2}=-\frac{\beta_{i-1}}{\Delta x}\left(\Phi_i-\Phi_{i-1}\right)
$$
which gives us for the central derivative in total:
$$
\frac{\partial}{\partial x}\left(D(x) \frac{\partial}{\partial x} \Phi(x)\right) \approx \frac{1}{\Delta x}\left(-B C_i \Phi_i-\frac{\beta_{i-1}}{\Delta x}\left(\Phi_i-\Phi_{i-1}\right)\right)
$$
therefore I need to set
$$
A_{1,1} = \frac{1}{dx}\left( -BC_{1}-\frac{D}{dx} \right)
$$
and
$$
A_{1,2} = - \frac{1}{dx^2}D
$$
# Implementing the Boundary Conditions
To implement the boundary conditions, we need to modify the first and last rows of the matrix \( M \) in order to enforce the specified boundary flux conditions at the left and right boundaries of the domain. Let's go through the details and code step-by-step.

### Problem Breakdown

The differential equation is
$$
D \frac{d^2}{dx^2} \phi - \Sigma_{a} \phi = 0,
$$
which is discretized using the finite difference method. The matrix \( M \) is constructed as:
$$
M = \text{streaming} + \text{collision},
$$
where:
- **streaming**: Represents the discretized Laplacian with the term \( D \frac{d^2}{dx^2} \).
- **collision**: Represents the reaction term \( - \Sigma_a \phi \).

### Boundary Conditions

1. **Left Boundary**:
   $$
   J^{-}(0) = \frac{\phi(0)}{4} + D \left. \frac{d\phi(x)}{dx} \right|_{x=0} = 0
   $$
   This implies:
   $$
   \frac{\phi(0)}{4} + D \frac{\phi(1) - \phi(0)}{dx} = 0.
   $$

2. **Right Boundary**:
   $$
   J^{+}(0) = \frac{\phi(0)}{4} - D \left. \frac{d\phi(x)}{dx} \right|_{x=0} = 0
   $$
   This implies:
   $$
   \frac{\phi(n-1)}{4} - D \frac{\phi(n-1) - \phi(n-2)}{dx} = 0.
   $$

### Modifying the Matrix \( M \)

To implement these boundary conditions, we adjust the first and last rows of the matrix \( M \) to incorporate the flux constraints. Let's break down the steps:

1. **Left Boundary**:
   - Modify the first row of \( M \) to enforce:
     $$
     \frac{\phi(0)}{4} + D \frac{\phi(1) - \phi(0)}{dx} = 0.
     $$
   - This means setting:
     $$
     M[1, 1] = \frac{1}{4} - \frac{D}{dx}, \quad M[1, 2] = \frac{D}{dx}.
     $$

2. **Right Boundary**:
   - Modify the last row of \( M \) to enforce:
     $$
     \frac{\phi(n-1)}{4} - D \frac{\phi(n-1) - \phi(n-2)}{dx} = 0.
     $$
   - This means setting:
     $$
     M[n, n] = \frac{1}{4} - \frac{D}{dx}, \quad M[n, n-1] = \frac{D}{dx}.
     $$

### Code Implementation

Here’s the Julia function to modify the matrix \( M \) based on the boundary conditions:

```julia
using SparseArrays

function apply_boundary_conditions!(M, D, dx, n)
    # Modify left boundary (J^-(0))
    M[1, 1] = 1/4 - D / dx
    M[1, 2] = D / dx

    # Modify right boundary (J^+(0))
    M[n, n] = 1/4 - D / dx
    M[n, n-1] = D / dx
end
```

### Usage

After constructing the initial matrix \( M \) as you showed, apply the boundary conditions by calling:

```julia
apply_boundary_conditions!(M, D, dx, n)
```

This function directly modifies the matrix \( M \) to reflect the boundary conditions on both ends of the domain.