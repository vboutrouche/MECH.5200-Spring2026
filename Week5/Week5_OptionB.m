%[text] # Activity: 2D Poisson Equation with Neumann Boundary Condition
%[text] **Objective:** Extend the 2D finite difference Poisson solver from Week 4 to handle a **Neumann** (flux) boundary condition on one edge, using the **ghost-node** technique.
%%
%[text] ## Background
%[text] In Week 4, we solved the 2D Poisson equation with **Dirichlet** (prescribed value) boundary conditions on all four edges. In many physical problems, we instead know the **flux** (derivative) at a boundary — this is a **Neumann** boundary condition.
%[text] ### Example: Steady-State Heat Conduction
%[text] Consider a rectangular plate $\[0,1\] \\times \[0,1\]$ with temperature $T(x,y)$ satisfying:
%[text]{"align":"center"} $\\frac{\\partial^2 T}{\\partial x^2} + \\frac{\\partial^2 T}{\\partial y^2} = -2\\pi^2 \\sin(\\pi x)\\sin(\\pi y)$
%[text] with boundary conditions:
%[text] - **Left** ($x=0$): $T = 0$ (Dirichlet)
%[text] - **Right** ($x=1$): $T = 0$ (Dirichlet)
%[text] - **Bottom** ($y=0$): $T = 0$ (Dirichlet)
%[text] - **Top** ($y=1$): $\\frac{\\partial T}{\\partial y} = -\\pi\\sin(\\pi x)$ **(Neumann — prescribed flux!)** \
%[text] The exact solution is $T(x,y) = \\sin(\\pi x)\\sin(\\pi y)$, enabling error verification.
%%
%[text] ## Ghost-Node Technique for Neumann BCs
%[text] At the top boundary ($j = M$), we want to impose $\\frac{\\partial T}{\\partial y}\\bigg|\_{i,M} = g\_i$.
%[text] **Step 1:** Imagine a "ghost" row of nodes at $j = M+1$ (outside the domain).
%[text] **Step 2:** Approximate the derivative with a central difference at $j = M$:
%[text]{"align":"center"} $\\frac{T\_{i,M+1} - T\_{i,M-1}}{2\\Delta y} = g\_i \\quad \\Longrightarrow \\quad T\_{i,M+1} = T\_{i,M-1} + 2\\Delta y \\cdot g\_i$
%[text] **Step 3:** Substitute into the standard 5-point stencil at $(i, M)$:
%[text]{"align":"center"} $\\frac{T\_{i-1,M} - 2T\_{i,M} + T\_{i+1,M}}{\\Delta x^2} + \\frac{T\_{i,M-1} - 2T\_{i,M} + T\_{i,M+1}}{\\Delta y^2} = f\_{i,M}$
%[text] Replacing $T\_{i,M+1}$:
%[text]{"align":"center"} $\\frac{T\_{i-1,M} - 2T\_{i,M} + T\_{i+1,M}}{\\Delta x^2} + \\frac{2T\_{i,M-1} - 2T\_{i,M}}{\\Delta y^2} = f\_{i,M} - \\frac{2g\_i}{\\Delta y}$
%[text] **Result:** The stencil at the Neumann boundary has:
%[text] - Diagonal: $-2/\\Delta x^2 - 2/\\Delta y^2$ (same as interior)
%[text] - West neighbor: $1/\\Delta x^2$ (same)
%[text] - East neighbor: $1/\\Delta x^2$ (same)
%[text] - South neighbor: $\\mathbf{2}/\\Delta y^2$ (**doubled!** because ghost node was eliminated)
%[text] - North neighbor: **gone** (no equation for ghost node)
%[text] - RHS: $f\_{i,M} - 2g\_i/\\Delta y$ \
%%
%[text] ## Step 1: Setup the Grid and Node Mapping
%[text] This is the same grid and mapping setup from Week 4.

N = 25;                    % Nodes in x-direction
M = 25;                    % Nodes in y-direction
x = linspace(0, 1, N);
y = linspace(0, 1, M);
dx = x(2) - x(1);
dy = y(2) - y(1);

% Node mapping: Map(ix, iy) gives the equation number for node (ix, iy)
Map = zeros(N, M);
Counter = 1;
for ix = 1:N
    for iy = 1:M
        Map(ix, iy) = Counter;
        Counter = Counter + 1;
    end
end
%%
%[text] ## Step 2: Define Source Term and Exact Solution
%[text] These functions will be used for assembling the RHS and for error checking.
%[text] - Source: $f(x,y) = -2\\pi^2 \\sin(\\pi x)\\sin(\\pi y)$
%[text] - Exact: $T(x,y) = \\sin(\\pi x)\\sin(\\pi y)$
%[text] - Neumann flux at top: $g(x) = \\frac{\\partial T}{\\partial y}\\bigg|\_{y=1} = -\\pi\\sin(\\pi x)$ \

source = @(xv, yv) -2*pi^2*sin(pi*xv).*sin(pi*yv);
T_exact = @(xv, yv) sin(pi*xv).*sin(pi*yv);
g_top = @(xv) -pi*sin(pi*xv);    % Neumann flux at y = 1
%%
%[text] ## Step 3: Assemble Interior Stencil
%[text] **TODO 1:** Loop over **interior** nodes $(ix = 2:N-1, \\; iy = 2:M-1)$ and fill the matrix with the standard 5-point Laplacian stencil. This is the same as Week 4.
%[text] **Reminder (5-point stencil):**
%[text]{"align":"center"} $A = \\pmatrix{ 0 & 1/\\Delta y^2 & 0 \\cr 1/\\Delta x^2 & -2/\\Delta x^2 - 2/\\Delta y^2 & 1/\\Delta x^2 \\cr 0 & 1/\\Delta y^2 & 0}$

NNodes = N*M;
A = zeros(NNodes);
b = zeros(NNodes, 1);

% --- YOUR CODE HERE (TODO 1) ---
% for ix = 2:N-1
%     for iy = 2:M-1
%         EqnId = Map(ix, iy);
%         A(EqnId, EqnId)          = ???   % Center
%         A(EqnId, Map(ix-1, iy))  = ???   % West
%         A(EqnId, Map(ix+1, iy))  = ???   % East
%         A(EqnId, Map(ix, iy-1))  = ???   % South
%         A(EqnId, Map(ix, iy+1))  = ???   % North
%         b(EqnId) = ???                    % Source term
%     end
% end
% --------------------------------
%%
%[text] ## Step 4: Apply Dirichlet Boundary Conditions
%[text] **TODO 2:** Apply Dirichlet BCs ($T = 0$) on the **left, right, and bottom** edges.
%[text] This is exactly the same pattern from Week 4: set the diagonal to 1 and the RHS to the boundary value.

% --- YOUR CODE HERE (TODO 2) ---
% Bottom (iy = 1): T = 0
% for ix = 1:N
%     ...
% end

% Left (ix = 1): T = 0
% for iy = 1:M
%     ...
% end

% Right (ix = N): T = 0
% for iy = 1:M
%     ...
% end
% --------------------------------
%%
%[text] ## Step 5: Apply Neumann BC at the Top
%[text] **TODO 3:** This is the **new** part! At the top boundary ($iy = M$), we impose $\\frac{\\partial T}{\\partial y} = g(x)$.
%[text] Use the ghost-node derivation from the Background section:
%[text] - For each interior column $ix = 2:N-1$ at $iy = M$:
%[text] -  - Center: $-2/\\Delta x^2 - 2/\\Delta y^2$
%[text] -  - West:   $1/\\Delta x^2$
%[text] -  - East:   $1/\\Delta x^2$
%[text] -  - South:  $2/\\Delta y^2$ (doubled!)
%[text] -  - RHS:    $f(x\_i, y\_M) - 2 g(x\_i) / \\Delta y$ \
%[text] **Note:** The corner nodes at $(1, M)$ and $(N, M)$ are already handled by the Dirichlet BCs on the left/right edges (from TODO 2).

% --- YOUR CODE HERE (TODO 3) ---
% for ix = 2:N-1
%     EqnId = Map(ix, M);
%     A(EqnId, EqnId)          = ???   % Center (same as interior)
%     A(EqnId, Map(ix-1, M))   = ???   % West
%     A(EqnId, Map(ix+1, M))   = ???   % East
%     A(EqnId, Map(ix, M-1))   = ???   % South (doubled!)
%     b(EqnId) = ???                    % Modified RHS
% end
% --------------------------------
%%
%[text] ## Step 6: Solve and Visualize
%[text] **TODO 4:** Solve the system $AT = b$, reshape the solution, and create a filled contour plot.

% --- YOUR CODE HERE (TODO 4) ---
% T = A\b;
% [X, Y] = meshgrid(x, y);
% T_grid = reshape(T, size(X));   % Hint: reshape to match meshgrid output
% contourf(X, Y, T_grid, 20, LineStyle="-")
% colorbar
% xlabel('x')
% ylabel('y')
% title('Temperature Distribution (Neumann BC at top)')
% --------------------------------
%%
%[text] ## Step 7: Error Analysis
%[text] **TODO 5:** Compute the exact solution on the grid and find the maximum error.

% --- YOUR CODE HERE (TODO 5) ---
% [X, Y] = meshgrid(x, y);
% T_ex = T_exact(X, Y);
% error = abs(T_grid - T_ex);
% fprintf('Max absolute error: %.2e\n', max(error(:)))
%
% contourf(X, Y, error, 20, LineStyle="-")
% colorbar
% xlabel('x')
% ylabel('y')
% title('Absolute Error')
% --------------------------------

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
