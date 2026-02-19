%[text] # Activity: 1D Finite Differences with Domain Mapping
%[text] **Objective:** Solve a 1D boundary value problem on a non-uniform physical grid by mapping it to a uniform computational domain.
%%
%[text] ## Background
%[text] So far, we have solved $\\frac{d^2 u}{dx^2} = f(x)$ on uniform grids. In many problems, the physical domain is not uniformly spaced, or we want to cluster nodes in regions of interest. **Domain mapping** allows us to solve on a uniform computational domain $\\xi \\in \[0,1\]$ while the physical grid $x(\\xi)$ can be non-uniform.
%[text] ### The Mapping
%[text] We define a mapping from computational coordinates $\\xi$ to physical coordinates $x$:
%[text]{"align":"center"} $x(\\xi) = -2\\cos(\\pi \\xi), \\quad \\xi \\in \[0,1\]$
%[text] This maps $\\xi = 0 \\to x = -2$ and $\\xi = 1 \\to x = 2$, with nodes clustered near the endpoints.
%[text] ### Transformed PDE
%[text] Using the chain rule, the second derivative in physical space becomes (in computational space):
%[text]{"align":"center"} $\\frac{d^2 u}{dx^2} = \\frac{1}{J^2}\\frac{d^2 u}{d\\xi^2} - \\frac{J'}{J^3}\\frac{du}{d\\xi}$
%[text] where $J = \\frac{dx}{d\\xi}$ is the **Jacobian** and $J' = \\frac{d^2 x}{d\\xi^2}$.
%[text] ### Mapped Finite Difference Stencil
%[text] Applying central differences in $\\xi$, the stencil for node $i$ becomes:
%[text]{"align":"center"} $A\_{i,i} = \\frac{-2}{J\_i^2 \\Delta\\xi^2}$
%[text]{"align":"center"} $A\_{i,i-1} = \\frac{1}{J\_i^2 \\Delta\\xi^2} + \\frac{J'\_i}{J\_i^3} \\cdot \\frac{1}{2\\Delta\\xi}$
%[text]{"align":"center"} $A\_{i,i+1} = \\frac{1}{J\_i^2 \\Delta\\xi^2} - \\frac{J'\_i}{J\_i^3} \\cdot \\frac{1}{2\\Delta\\xi}$
%%
%[text] ## Problem Statement
%[text] Solve the following 1D BVP using domain mapping:
%[text]{"align":"center"} $\\frac{d^2 u}{dx^2} = -\\pi^2 \\sin(\\pi x)$
%[text] on $x \\in \[-2, 2\]$ with Dirichlet boundary conditions:
%[text]{"align":"center"} $u(-2) = 0, \\quad u(2) = 0$
%[text] The exact solution is $u(x) = \\sin(\\pi x)$, which you will use to verify your code.
%%
%[text] ## Step 1: Setup the Computational Grid
%[text] We work on a uniform grid in $\\xi \\in \[0,1\]$.

N = 51;                % Number of nodes
dxi = 1/(N-1);         % Uniform spacing in computational domain
xi = 0:dxi:1;          % Computational grid
%%
%[text] ## Step 2: Compute the Mapping
%[text] **TODO 1:** Compute the physical coordinates $x(\\xi)$, the Jacobian $J = \\frac{dx}{d\\xi}$, and the second derivative of the mapping $J' = \\frac{d^2x}{d\\xi^2}$.
%[text] Recall:
%[text] - $x(\\xi) = -2\\cos(\\pi\\xi)$
%[text] - $J(\\xi) = \\frac{dx}{d\\xi} = 2\\pi\\sin(\\pi\\xi)$
%[text] - $J'(\\xi) = \\frac{d^2x}{d\\xi^2} = 2\\pi^2\\cos(\\pi\\xi)$ \

% --- YOUR CODE HERE (TODO 1) ---
% x = ???
% J = ???
% Jp = ???   (J' = d2x/dxi2)
% --------------------------------
%%
%[text] ## Step 3: Assemble the Matrix
%[text] **TODO 2:** Loop over interior nodes ($i = 2$ to $N-1$) and fill the sparse matrix $A$ and the right-hand side vector $f$ using the mapped stencil formulas from the Background section.
%[text] **Hint:** The RHS at each node is $f\_i = -\\pi^2 \\sin(\\pi x\_i)$, evaluated at the **physical** coordinate.

A = spalloc(N, N, 3*N);   % Sparse tridiagonal matrix
f = zeros(N, 1);           % RHS vector

% --- YOUR CODE HERE (TODO 2) ---
% for i = 2:(N-1)
%     A(i,i)   = ???
%     A(i,i-1) = ???
%     A(i,i+1) = ???
%     f(i)     = ???
% end
% --------------------------------
%%
%[text] ## Step 4: Apply Boundary Conditions
%[text] **TODO 3:** Apply Dirichlet boundary conditions: $u(-2) = 0$ and $u(2) = 0$.
%[text] This is the same pattern you used in previous activities: set the diagonal to 1 and the RHS to the boundary value.

% --- YOUR CODE HERE (TODO 3) ---
% Left boundary (xi = 0, x = -2): u = 0
% Right boundary (xi = 1, x = 2): u = 0
% --------------------------------
%%
%[text] ## Step 5: Solve and Visualize
%[text] **TODO 4:** Solve the system $Au = f$ and create a 2x2 subplot showing:
%[text] - (1,1): Computational domain nodes ($\\xi$ vs 0)
%[text] - (1,2): Physical domain nodes ($x$ vs 0)
%[text] - (2,1): Solution in computational domain ($\\xi$ vs $u$)
%[text] - (2,2): Solution in physical domain ($x$ vs $u$) \

% --- YOUR CODE HERE (TODO 4) ---
% u = ???
% subplot(2,2,1)
%   plot the computational grid
% subplot(2,2,2)
%   plot the physical grid
% subplot(2,2,3)
%   plot u vs xi
% subplot(2,2,4)
%   plot u vs x
% --------------------------------
%%
%[text] ## Step 6: Error Analysis
%[text] **TODO 5:** Compute the exact solution at the physical grid points and plot the error.
%[text] The exact solution is $u\_{exact}(x) = \\sin(\\pi x)$.

% --- YOUR CODE HERE (TODO 5) ---
% u_exact = ???
% error = ???
% plot(x, error, 'r-o', 'LineWidth', 1.5)
% xlabel('x')
% ylabel('Error')
% title('Error: Numerical - Exact')
% grid on
% --------------------------------

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
