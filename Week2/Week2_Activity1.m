%[text] # Activity 1: Finite Difference Approximations
%[text] **Learning Objectives:** Implement and compare forward, backward, and central difference formulas for approximating derivatives.
%%
%[text] ## Step 1 - Set the problem
%[text] ### Setup
%[text] Approximate the derivative of one of the following using different methods:
%[text] - $ f(x) = \\sin(2\\pi x) $ at $ x = 0.1 $
%[text] - $f(x) = x^3$ at $x = 0$
%[text] - $f(x) = \\exp(x)$ at $x= 1$
%[text] - $f(x) = \\ln(x)$ at $x = 3$
%[text] - $f(x) = \\frac{1}{x}$ at $x=\\frac{1}{2}$
%[text] - $f(x) = \\tan(x)$ at $x=\\frac{1}{2}$
%[text] - $f(x) = x\\exp(x)$ at $x = -0.1$ \
%%
%[text] ### Forward Difference
%[text] The forward difference uses current and next point: 
%[text]{"align":"center"} $\\frac{\\partial u}{\\partial x}\\bigg|i = \\frac{u\_{i+1} - u\_i}{\\Delta x} + O(\\Delta x)$
%%
%[text] ### Backward Difference
%[text] The backward difference uses current and previous point: 
%[text]{"align":"center"} $\\frac{\\partial u}{\\partial x}\\bigg|i = \\frac{u\_i - u\_{i-1}}{\\Delta x} + O(\\Delta x)$ 
%%
%[text] ### Central Difference
%[text] The central difference uses both sides: 
%[text]{"align":"center"} $\\frac{\\partial u}{\\partial x}\\bigg|i = \\frac{u\_{i+1} - u\_{i-1}}{2\\Delta x} + O(\\Delta x^2)$ 
%%
%[text] ### Comparison
%[text] Compare the numerical approximation to the analytical solution.
%%
%[text] ## Step 2 - Convergence study
%[text] Modify the code above (or below) to check perform a convergence study. In other words, study how the error increases/reduces as the $\\Delta x$ is changing.
%%
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
