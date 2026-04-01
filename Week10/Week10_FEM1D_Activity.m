%[text] # Week 10 Solution: Building a 1D FEM Solver for String Deflection
%[text] This is the **instructor solution** file. It is structurally identical to the Activity file. Every TODO block is replaced with the correct code. Use this to check student work or project the answer after each TODO is attempted.
%%
%[text] ## 0) Problem Statement and Physical Setup
%[text] We solve the steady-state string deflection problem:
%[text]{"align":"center"} $ -u''(x) = f(x),\\quad x\\in(0,1)\\\\ u(0)=u(1)=0 $
%[text] Here $ u(x) $ is vertical deflection and $ f(x) $ is the transverse distributed load. We will use uniform load $ f(x)=1 $ and domain length $ L=1 $.
%[text] **FEM roadmap:** (1) derive the weak form, (2) define basis functions on a mesh, (3) compute the element stiffness matrix, (4) compute the element load vector, (5) assemble the global system, (6) apply BCs and solve — each step follows directly from the one before.

L        = 1.0;
f_handle = @(x) 1 + 0*x;

x_plot = linspace(0, L, 300);

figure
plot(x_plot, f_handle(x_plot), 'LineWidth', 1.8);
grid on;
xlabel('x');
ylabel('f(x)');
title('Uniform Distributed Load  f(x) = 1');
%%
%[text] ## 1) Weak Form via Integration by Parts
%[text] Direct discretization of $ -u''=f $ requires second-order finite-difference stencils. FEM avoids this: we multiply both sides by a smooth test function $ v $ satisfying $ v(0)=v(1)=0 $, integrate over the domain, then use integration by parts to transfer one derivative from $ u $ onto $ v $. The result involves only **first derivatives** — which piecewise-linear polynomials can provide.
%[text] Multiply both sides of $ -u''=f $ by a test function $ v $ satisfying $ v(0)=v(1)=0 $, and integrate:
%[text]{"align":"center"} $ \\int\_0^1 (-u'')v,dx = \\int\_0^1 fv,dx $
%[text] Apply the product rule $ \\frac{d}{dx}(u'v) = u''v + u'v' $ to trade one derivative from $ u $ onto $ v $:
%[text]{"align":"center"} $ \\int\_0^1 u'v',dx - \\left\[u'v \\right\]\_0^1 = \\int\_0^1 fv,dx $
%[text] Because $ v(0)=v(1)=0 $ the boundary term vanishes. The **weak form** is:
%[text]{"align":"center"} $ \\int\_0^1 u'v',dx = \\int\_0^1 fv,dx $
%[text] This replaces a second-order ODE with an integral equation that only needs first derivatives. We now approximate the infinite-dimensional solution $ u $ with piecewise-linear polynomials defined on a finite mesh — this leads to the finite linear system $ KU=F $ that we will assemble in sections 3–5.

%%
%[text] ## 2) Linear Basis Functions on the Reference Element
%[text] Partition $ \[0,L\] $ into $ n\_{el} $ elements. Map each physical element $ \[x\_a,x\_b\] $ to the reference domain $ \\xi\\in\[-1,1\] $ via:
%[text] $ x(\\xi)=\\frac{x\_a+x\_b}{2}+\\frac{h\_e}{2}\\xi,\\quad J=\\frac{dx}{d\\xi}=\\frac{h\_e}{2} $
%[text] The two linear shape functions on the reference element are:
%[text] $ N\_1(\\xi)=\\frac{1-\\xi}{2},\\quad N\_2(\\xi)=\\frac{1+\\xi}{2} $
%[text] Their $ \\xi $-derivatives are constant, meaning the approximated strain is uniform within each element.

syms xi_s real
N1s = (1 - xi_s)/2;
N2s = (1 + xi_s)/2;

xi_num = linspace(-1, 1, 200);
N1_num = double(subs(N1s, xi_s, xi_num));
N2_num = double(subs(N2s, xi_s, xi_num));

figure
plot(xi_num, N1_num, 'b-', 'LineWidth', 1.8);
hold on;
plot(xi_num, N2_num, 'r-', 'LineWidth', 1.8);

grid on;
xlabel('\xi');
ylabel('Value');
title('Reference Basis Functions and Their \xi-Derivatives');
legend('N_1', 'N_2', 'Location', 'bestoutside');
%%
%[text] ## 3) Element Stiffness Matrix — Symbolic Integration
%[text] On one element with nodes at $ x\_a $ and $ x\_b $, we approximate $ u\_h = N\_1 U\_1 + N\_2 U\_2 $. Substituting into the weak form and choosing $ v = N\_a $ for each row — the **Galerkin choice**, using the same shape functions as both trial and test functions — the element weak form becomes a 2×2 linear system:
%[text] $ k\_{ab}^{(e)}\\,U\_b = f\_a^{(e)}, \\quad\\text{where}\\quad k\_{ab}^{(e)} = \\int\_{x\_a}^{x\_b} N\_a'(x)\\,N\_b'(x),dx $
%[text] Mapping to the reference element $ \\xi\\in\[-1,1\] $ with chain rule $ dN\_a/dx = (2/h\_e)\\,dN\_a/d\\xi $ and Jacobian $ J = h\_e/2 $, this becomes:
%[text]{"align":"center"} $ k\_{ab}^{(e)} = \\int\_{-1}^{1} \\frac{dN\_a}{dx}\\frac{dN\_b}{dx}\\,J\\,d\\xi $
%[text] Evaluating all four entries ($ a,b\\in\\{1,2\\} $) gives the 2×2 element stiffness matrix.


%%
%[text] ## 4) Element Load Vector — Symbolic Integration
%[text] The element system $ k^{(e)} U^{(e)} = f^{(e)} $ has two parts. Section 3 built the **left-hand side** stiffness $ k^{(e)} $. Now we build the **right-hand side** $ f^{(e)} $: it distributes the applied load $ f(x) $ to the two element nodes by integrating $ f $ against each shape function.
%[text] Map the forcing to the reference coordinate: $ x(\\xi)=\\frac{x\_1+x\_2}{2}+\\frac{h\_e}{2}\\xi $. For affine forcing $ f=f\_0$, the load integral for node $ a $ is:
%[text] $ f\_a^{(e)} = \\int\_{-1}^{1} f(x(\\xi))\\,N\_a(\\xi)\\,J\\,d\\xi,\\quad a=1,2 $

%%
%[text] ## 5) Assemble the Global System
%[text] Sections 3 and 4 derived $ k^{(e)} $ and $ f^{(e)} $ for a **single** element. Now we repeat this for every element on the mesh and accumulate all contributions into the global $ n\\times n $ system $ KU=F $.
%[text] The rule for element $ e $ with global node indices `nodes = [i, j]`: scatter `ke` into `K(nodes,nodes)` and `fe` into `F(nodes)`. These two lines of MATLAB — repeated for each element — encode all the mesh connectivity.

%[text] The assembled $ K $ is symmetric and tridiagonal, but it is **singular** — no displacement is fixed yet, so the system has infinitely many solutions. Section 6 imposes the Dirichlet boundary conditions to make the system uniquely solvable.
%%
%[text] ## 6) Apply Boundary Conditions and Solve
%[text] The assembled $ K $ is singular — it needs boundary conditions to have a unique solution. The simplest approach is **row replacement**: wipe out the assembled equation for each boundary node and substitute the constraint directly.
%[text] For node 1 (at $ x=0 $): zero out row 1, set $ K(1,1)=1 $, set $ F(1)=0 $. Row 1 now reads $ 1 \\cdot U\_1 = 0 $, i.e. $ u(0)=0 $. Repeat for the last node: $ K(n,:)=0 $, $ K(n,n)=1 $, $ F(n)=0 $ → $ u(1)=0 $.
%[text] The modified $ K $ is no longer singular. Solving $ U = K \\backslash F $ gives the complete nodal solution with both BCs enforced.

%%
%[text] ## 7) Validation and Convergence Study
%[text] With the system solved, we verify accuracy in two ways: (1) compare the nodal solution against the known exact answer for $ f=1 $, and (2) run a mesh-refinement study to confirm the expected $ h^2 $ convergence rate.
%[text] For $ f=1 $, $ L=1 $, the exact solution is $ u=x(1-x)/2 $. Because this is a quadratic and the Gauss quadrature is exact, nodal FEM values are recovered to machine precision. The convergence plot therefore uses $ f=\\pi^2\\sin(\\pi x) $ with exact $ u=\\sin(\\pi x) $, where the FEM has genuine $ h^2 $ error.




%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
