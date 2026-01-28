%[text] # Activity 2: Second Derivatives and String Problem Setup
%[text] Learn the second derivative finite difference formula and apply it to the string deflection boundary value problem.
%%
%[text] ## Second Derivative Formula
%[text] Approximate the second derivative using three neighboring points: 
%[text]{"align":"center"} $\\frac{\\partial^2 u}{\\partial x^2}\\bigg|i = \\frac{u{i+1} - 2u\_i + u\_{i-1}}{\\Delta x^2} + O(\\Delta x^2)$ 
%[text] This is second-order accurate.

%%
%[text] ## Physical Problem: String Deflection BVP
%[text] A string under distributed load satisfies: $ \\frac{d^2 u}{dx^2} = f(x) $ on $ \[0,1\] $ with $ u(0) = 0 $ and $ u(1) = 0 $ (fixed ends). Here $ u(x) $ is deflection and $ f(x) $ is the load.


%%
%[text] ## Discretized Linear System
%[text] At each interior point, the PDE becomes: 
%[text]{"align":"center"} $ u\_{i+1} - 2u\_i + u\_{i-1} = f\_i (\\Delta x)^2 $. 
%[text] This gives $ (n-2) $ equations in $ (n-2) $ unknowns with a tridiagonal matrix structure.



%[text] **Key insight:** The tridiagonal structure (sparse matrix) makes this system very efficient to solve, even for millions of grid points.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
