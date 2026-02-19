%[text] # Activity: 2D Bilinear Mapping & Visualization
%[text] **Objective:** Map a unit computational square $(\\xi, \\eta) \\in \[0,1\]^2$ to an irregular physical quadrilateral using a **bilinear mapping**, visualize the transformation, and analyze the Jacobian.
%%
%[text] ## Background
%[text] In many engineering problems, the physical domain is not a simple rectangle. A **bilinear mapping** allows us to transform a unit square in computational space to a general quadrilateral (e.g., a trapezoidal duct cross-section) in physical space.
%[text] ### Bilinear Mapping Formula
%[text] The mapping from $(\\xi, \\eta)$ to $(x, y)$ is:
%[text]{"align":"center"} $x(\\xi, \\eta) = R\_x \\xi\\eta + S\_x \\xi + T\_x \\eta + U\_x$
%[text]{"align":"center"} $y(\\xi, \\eta) = R\_y \\xi\\eta + S\_y \\xi + T\_y \\eta + U\_y$
%[text] The four coefficients $(R, S, T, U)$ for each coordinate are determined by evaluating the mapping at the four corners of the unit square:
%[text:table]
%[text] | Corner | $\\xi$ | $\\eta$ | Physical $(x, y)$ |
%[text] | --- | --- | --- | --- |
%[text] | 1 | 0 | 0 | $(x\_1, y\_1)$ |
%[text] | 2 | 1 | 0 | $(x\_2, y\_2)$ |
%[text] | 3 | 1 | 1 | $(x\_3, y\_3)$ |
%[text] | 4 | 0 | 1 | $(x\_4, y\_4)$ |
%[text:table]
%[text] 
%[text] Substituting each corner into the formula gives a $4 \\times 4$ system:
%[text]{"align":"center"} $\\pmatrix{ \\xi\_1\\eta\_1 & \\xi\_1 & \\eta\_1 & 1 \\cr \\xi\_2\\eta\_2 & \\xi\_2 & \\eta\_2 & 1 \\cr \\xi\_3\\eta\_3 & \\xi\_3 & \\eta\_3 & 1 \\cr \\xi\_4\\eta\_4 & \\xi\_4 & \\eta\_4 & 1 } \\pmatrix{ R \\cr S \\cr T \\cr U } = \\pmatrix{ x\_1 \\cr x\_2 \\cr x\_3 \\cr x\_4 }$
%[text] Which simplifies to:
%[text]{"align":"center"} $\\pmatrix{ 0 & 0 & 0 & 1 \\cr 0 & 1 & 0 & 1 \\cr 1 & 1 & 1 & 1 \\cr 0 & 0 & 1 & 1 } \\pmatrix{ R \\cr S \\cr T \\cr U } = \\pmatrix{ x\_1 \\cr x\_2 \\cr x\_3 \\cr x\_4 }$
%[text] ### Jacobian of the Mapping
%[text] The Jacobian matrix describes how the mapping stretches/compresses space:
%[text]{"align":"center"} $\\mathbf{J} = \\pmatrix{ \\frac{\\partial x}{\\partial \\xi} & \\frac{\\partial x}{\\partial \\eta} \\cr \\frac{\\partial y}{\\partial \\xi} & \\frac{\\partial y}{\\partial \\eta} }$
%[text] The **determinant** $|\\mathbf{J}|$ measures local area scaling. If $|\\mathbf{J}| \> 0$ everywhere, the mapping is valid (no folding).
%%
%[text] ## Physical Context
%[text] We'll map the unit square to a **trapezoidal duct cross-section** with corners:
%[text] - Corner 1: $(0, 0)$
%[text] - Corner 2: $(2, 0)$
%[text] - Corner 3: $(1.5, 1)$
%[text] - Corner 4: $(0.5, 1)$ \
%%
%[text] ## Step 1: Define the Physical Domain Corners

% Physical domain corner coordinates (trapezoidal cross-section)
% Corner ordering: (0,0) -> (1,0) -> (1,1) -> (0,1) in (xi,eta)
corners_x = [0.0; 2.0; 1.5; 0.5];   % x-coordinates of 4 corners
corners_y = [0.0; 0.0; 1.0; 1.0];   % y-coordinates of 4 corners
%%
%[text] ## Step 2: Build and Solve the Mapping System
%[text] **TODO 1:** Set up the $4 \\times 4$ coefficient matrix and solve for the mapping coefficients $\[R, S, T, U\]$ for both $x$ and $y$.
%[text] **Hint:** Each row of the matrix corresponds to one corner evaluated at $(\\xi\_i, \\eta\_i)$. The columns represent $\[\\xi\\eta, \\; \\xi, \\; \\eta, \\; 1\]$.

% --- YOUR CODE HERE (TODO 1) ---
% Coeff = ???    % 4x4 matrix (see formula above)
% RSTU_x = ???   % Solve for x mapping coefficients
% RSTU_y = ???   % Solve for y mapping coefficients
% --------------------------------

% Display the coefficients
% fprintf('x-mapping: R=%.3f, S=%.3f, T=%.3f, U=%.3f\n', RSTU_x)
% fprintf('y-mapping: R=%.3f, S=%.3f, T=%.3f, U=%.3f\n', RSTU_y)
%%
%[text] ## Step 3: Create the Mapped Grid
%[text] **TODO 2:** Create a uniform grid in $(\\xi, \\eta)$ using `meshgrid`, then compute the corresponding physical coordinates $(x, y)$ using the bilinear mapping formula.

N = 21;                          % Grid resolution
dh = 1/(N-1);
xi  = 0:dh:1;
eta = 0:dh:1;
[xiM, etaM] = meshgrid(xi, eta);

% --- YOUR CODE HERE (TODO 2) ---
% xM = ???   % Physical x-coordinates (use RSTU_x and the bilinear formula)
% yM = ???   % Physical y-coordinates (use RSTU_y and the bilinear formula)
% --------------------------------
%%
%[text] ## Step 4: Plot Both Domains Side-by-Side
%[text] **TODO 3:** Create a side-by-side comparison of the computational and physical grids.

% --- YOUR CODE HERE (TODO 3) ---
% subplot(1,2,1)
%   Plot the computational grid (xi vs eta)
%   - Draw the unit square boundary: plot([0 1 1 0 0], [0 0 1 1 0], ...)
%   - Plot all grid points: plot(xiM, etaM, ...)
%   - Label axes, add title
%
% subplot(1,2,2)
%   Plot the physical grid (xM vs yM)
%   - Draw the quadrilateral boundary:
%     plot([corners_x; corners_x(1)], [corners_y; corners_y(1)], ...)
%   - Plot all grid points: plot(xM, yM, ...)
%   - Label axes, add title
%   - axis equal
% --------------------------------
%%
%[text] ## Step 5: Compute and Visualize the Jacobian
%[text] **TODO 4:** Compute the four partial derivatives at each grid point.
%[text] From the bilinear mapping formula:
%[text] - $\\frac{\\partial x}{\\partial \\xi} = R\_x \\eta + S\_x$
%[text] - $\\frac{\\partial x}{\\partial \\eta} = R\_x \\xi + T\_x$
%[text] - $\\frac{\\partial y}{\\partial \\xi} = R\_y \\eta + S\_y$
%[text] - $\\frac{\\partial y}{\\partial \\eta} = R\_y \\xi + T\_y$ \

% --- YOUR CODE HERE (TODO 4) ---
% dxdxi  = ???
% dxdeta = ???
% dydxi  = ???
% dydeta = ???
% --------------------------------
%%
%[text] ## Step 6: Compute and Plot the Jacobian Determinant
%[text] **TODO 5:** Compute $|\\mathbf{J}| = \\frac{\\partial x}{\\partial \\xi}\\frac{\\partial y}{\\partial \\eta} - \\frac{\\partial x}{\\partial \\eta}\\frac{\\partial y}{\\partial \\xi}$ and plot it on the physical domain.
%[text] **Interpretation:** $|\\mathbf{J}|$ tells you how much the mapping stretches or compresses local area. If $|\\mathbf{J}|$ varies significantly, the mesh quality is non-uniform.

% --- YOUR CODE HERE (TODO 5) ---
% detJ = ???
%
% contourf(xM, yM, detJ, 20, LineStyle="-")
% colorbar
% axis equal
% xlabel('x')
% ylabel('y')
% title('Jacobian Determinant |J| on Physical Domain')
%
% fprintf('Min |J| = %.4f, Max |J| = %.4f\n', min(detJ(:)), max(detJ(:)))
% --------------------------------
%%
%[text] ## Bonus: Solve the Poisson Equation on the Mapped Domain
%[text] This section is **optional** for those who finish early. Solve:
%[text]{"align":"center"} $\\nabla^2 u = 1$ (membrane deflection)
%[text] on the trapezoidal domain with $u = 0$ on all boundaries (Dirichlet).
%[text] **The transformed Laplacian in computational coordinates** (for a general mapping) is:
%[text]{"align":"center"} $\\nabla^2 u = \\frac{1}{|J|}\\left\[\\frac{\\partial}{\\partial \\xi}\\left(\\frac{\\alpha}{|J|}\\frac{\\partial u}{\\partial \\xi}\\right) - \\frac{\\partial}{\\partial \\xi}\\left(\\frac{\\beta}{|J|}\\frac{\\partial u}{\\partial \\eta}\\right) - \\frac{\\partial}{\\partial \\eta}\\left(\\frac{\\beta}{|J|}\\frac{\\partial u}{\\partial \\xi}\\right) + \\frac{\\partial}{\\partial \\eta}\\left(\\frac{\\gamma}{|J|}\\frac{\\partial u}{\\partial \\eta}\\right)\\right\]$
%[text] where:
%[text] - $\\alpha = \\left(\\frac{\\partial x}{\\partial \\eta}\\right)^2 + \\left(\\frac{\\partial y}{\\partial \\eta}\\right)^2$
%[text] - $\\beta = \\frac{\\partial x}{\\partial \\xi}\\frac{\\partial x}{\\partial \\eta} + \\frac{\\partial y}{\\partial \\xi}\\frac{\\partial y}{\\partial \\eta}$
%[text] - $\\gamma = \\left(\\frac{\\partial x}{\\partial \\xi}\\right)^2 + \\left(\\frac{\\partial y}{\\partial \\xi}\\right)^2$ \
%[text] **Hint:** For the trapezoidal domain with the given corners, $\\beta = 0$ when the mapping is orthogonal. Check the value of $\\beta$ at a few points — if it's small, you can simplify by ignoring the cross-derivative terms.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
