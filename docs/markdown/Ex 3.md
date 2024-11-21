

$$
\frac{1}{dx^2}(\phi_{1} - 2\phi_{2} + \phi_{3}) = C
$$
$$
\frac{1}{dx^2} \begin{pmatrix}
-2  & 1 \\
1 & -2 & 1 \\
 & 1 & -2 
\end{pmatrix} \begin{pmatrix}
\phi_{2} \\
\phi_{3} \\
\phi_{4}
\end{pmatrix} = \begin{pmatrix}
c_{1} - \frac{\phi_{1}}{dx^2}\\
c_{2} \\
c_{3} - \frac{\phi_{5}}{dx^2}
\end{pmatrix}
$$
where $\phi_{1}$ and $\phi_{5}$ are the dirichlet boundary conditions of the problem

# Question 1
%% Relationship between $\mathrm{\Phi}_i,\ \mathrm{\Phi}_{i+1},\ \mathrm{\Phi}_{i-1}$ for group g at any point within the material: %%
Using the laplace operator, that was implemented in the previous exercises
$$
\triangle \phi_{i} = \frac{1}{dx^2}(\phi_{i+1}- 2\phi_{i} + \phi_{i-1})
$$
the diffusion equation with downscattering and fission can be expressed at any point inside of the material as
$$
\begin{align}
D_{1}\triangle\phi_{1}-\Sigma_{a}\phi_{1} - \Sigma_{1\to2}\phi_{1}  & = -\frac{1}{k}(\nu \Sigma_{f_{1}}\phi_{1} + \nu \Sigma_{f_{2}}\phi_{2}) \\
D_{2} \triangle\phi_{2} - \Sigma_{a}\phi_{2} + \Sigma_{1\to_{2}} \phi_{2}  & = 0
\end{align}
$$
where $\phi_{1}$ and $\phi_{2}$ are understood to be at point $i$.
The two equations can be organized in this matrix.
$$
\begin{pmatrix}
D_{s} \triangle - \Sigma_{a} - \Sigma_{1\to2}  & 0 \\
\Sigma_{1\to2}  & D_{f} \triangle - \Sigma_{a}
\end{pmatrix} \begin{pmatrix}
\phi_{f} \\
\phi_{s}
\end{pmatrix} = -
\frac{1}{k} \begin{pmatrix}
\nu\Sigma_{ff}  & \nu\Sigma_{fs} \\
0 & 0
\end{pmatrix} \begin{pmatrix}
\phi_{f} \\
\phi_{s}
\end{pmatrix}
$$
# Question 2 - Analytical Solution
We make a seperation of variables, seperating $\phi_{g}(x)$ into an energy independant $\phi(x)$ and an energy group $g$ dependant amplitude $\varphi_{g}$:
$$
\phi_{g}(x) = \varphi_{g}\phi(x)
$$
then using that
$$
\triangle \phi(x) = B\phi(x)
$$
we can rewrite the equation
$$
\begin{pmatrix}
D_{s} B - \Sigma_{a} - \Sigma_{1\to2}  & 0 \\
\Sigma_{1\to2}  & D_{f} B - \Sigma_{a}
\end{pmatrix} \begin{pmatrix}
\varphi_{f} \\
\varphi_{s}
\end{pmatrix} = -
\frac{1}{k} \begin{pmatrix}
\nu\Sigma_{ff}  & \nu\Sigma_{fs} \\
0 & 0
\end{pmatrix} \begin{pmatrix}
\varphi_{f} \\
\varphi_{s}
\end{pmatrix}
$$
which is $A\varphi = \frac{F}{k}\varphi$
since A is a 2x2 matrix it can be easily inverted:
$$
A^{-1} = \frac{1}{\det(A)} \begin{pmatrix} (D_{f} B - \Sigma_{a}) & 0 \\ -\Sigma_{1\to2} & D_{s} B - \Sigma_{a} - \Sigma_{1\to2} \end{pmatrix}
$$
The determinant of \( A \) is given by:
$$
\det(A) = (D_{s} B - \Sigma_{a} - \Sigma_{1\to2})(D_{f} B - \Sigma_{a}) 
$$
