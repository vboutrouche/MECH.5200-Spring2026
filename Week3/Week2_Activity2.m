%[text] # Activity 2: Second Derivatives and String Problem Setup
%[text] Learn the second derivative finite difference formula and apply it to the string deflection boundary value problem.
%%
%[text] ## Second Derivative Formula
%[text] Approximate the second derivative using three neighboring points: 
%[text]{"align":"center"} $\\frac{\\partial^2 u}{\\partial x^2}\\bigg|i = \\frac{u{i+1} - 2u\_i + u\_{i-1}}{\\Delta x^2} + O(\\Delta x^2)$ 
%[text] This is second-order accurate.
%%
%[text] ## Physical Problem: String Deflection BVP
%[text] A string under distributed load satisfies: $ \\frac{d^2 u}{dx^2} = 1 $ on $ \[0,1\] $ with $ u(0) = 0 $ and $ u(1) = 0 $ (fixed ends). Here $ u(x) $ is deflection and we have a constant unit load $ f(x) = 1 $.
%%
%[text] ## Discretized Linear System
%[text] At each interior point, the PDE becomes: 
%[text]{"align":"center"} $ u\_{i+1} - 2u\_i + u\_{i-1} = (\\Delta x)^2 $. 
%[text] This gives $ (n-2) $ equations in $ (n-2) $ unknowns with a tridiagonal matrix structure.
%%
%[text] ## Matrix-Vector Format
%[text] Rearranging the discretized equations: $ u\_{i+1} - 2u\_i + u\_{i-1} = (\\Delta x)^2 $ for $ i = 1, 2, \\ldots, n-2 $
%[text] We can write this system in matrix form: 
%[text]{"align":"center"} $ \\mathbf{A} \\mathbf{u} = \\mathbf{b} $
%[text] where $ \\mathbf{u} = \[u\_1, u\_2, \\ldots, u\_{n-2}\]^T $ contains the unknowns at interior points, and $ \\mathbf{b} = \[(\\Delta x)^2, (\\Delta x)^2, \\ldots, (\\Delta x)^2\]^T $.
%[text] 
%[text] The coefficient matrix $ \\mathbf{A} $ is **tridiagonal**:
%[text]{"align":"center"} $ \\mathbf{A} = \\pmatrix{ -2 & 1 & 0 & \\cdots & 0 \\cr 1 & -2 & 1 & \\cdots & 0 \\cr 0 & 1 & -2 & \\cdots & 0 \\cr \\vdots & \\vdots & \\vdots & \\ddots & \\vdots \\cr 0 & 0 & 0 & \\cdots & -2 }$


%[text] **Key insight:** The tridiagonal structure (sparse matrix) makes this system very efficient to solve, even for millions of grid points.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
