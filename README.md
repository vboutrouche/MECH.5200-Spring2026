# MECH.5200 - Numerical Methods in Computational Mechanics

Welcome to the course materials for MECH.5200. This repository contains teaching scripts, code examples, and activities for learning numerical methods in computational mechanics.

## Overview

This course covers fundamental numerical techniques used in computational mechanics, with emphasis on practical implementation and understanding of algorithms through hands-on coding.

## Repository Structure

### Week 1: Foundations
**Location:** `Week1/`

Introduction to numerical computing fundamentals:
- **Vectors and Linear Algebra Operations**
  - `multiplicationInner.mlx` - Inner product computations
  - `multiplactionVector.mlx` - Vector multiplication operations

- **Introductory Scripts**
  - `script1.m` - Basic MATLAB programming concepts
  - `script2.m` - Mathematical operations fundamentals
  - `script3.m` - Advanced vector/matrix operations

- **Slides:** Lecture presentation materials

### Week 2: Polynomial Approximation Methods
**Location:** `Week2/`

Explore different approaches to function approximation and interpolation:

- **`Lagrange.m`** - Lagrange Polynomial Interpolation
  - Demonstrates how to construct polynomials that pass exactly through given data points
  - Uses Lagrange basis functions: $L_j(x_j) = 1$ and $L_j(x_m) = 0$ for $m \neq j$
  - Ideal for understanding interpolation theory and implementation
  - **Key Concept:** Direct polynomial construction without solving linear systems

- **`TaylorExpansion.m`** - Taylor Series Approximation  
  - Shows how to approximate functions using Taylor polynomials
  - Compares original functions with truncated series expansions
  - Demonstrates convergence of higher-order terms
  - **Key Concept:** Local approximation based on derivatives at an expansion point

- **Activities**
  - `Week2_Activity1.m` - Practical exercises on interpolation
  - `Week2_Activity2.m` - Applications and comparisons of approximation methods

## Script Format

All teaching scripts follow the **MATLAB Live Script** format (plain text `.m` files with rich formatting):

- **Rich Text Sections:** Use `%[text]` markup for explanations
- **Mathematical Equations:** LaTeX formatting with inline ($...$) and display math support
- **Code Structure:** 
  - `%%` markers separate logical sections
  - Clear variable naming aligned with mathematical notation
  - Detailed comments explaining algorithms
- **Visualization:** Implicit figure creation with `fplot()` and `scatter()` commands

## Learning Objectives

By working through these materials, you will understand:

1. **Polynomial Interpolation**
   - Lagrange basis function construction
   - Error analysis and convergence properties
   - Practical applications in numerical analysis

2. **Function Approximation**
   - Taylor series theory and implementation
   - Truncation error and remainder terms
   - Convergence properties of series expansions

3. **Numerical Implementation**
   - Efficient MATLAB coding for scientific computing
   - Visualization of mathematical concepts
   - Comparison of different approximation methods

## Using These Materials

### Running a Script
```matlab
% Open the file in MATLAB
open('Lagrange.m')

% Run the entire script or individual sections
```

### Modifying and Experimenting
- Change `N` (number of points) to see interpolation quality change
- Adjust the `Interval` to explore different domains
- Try different functions (replace `sin(x)*exp(x)` with other expressions)

### Best Practices
- Run scripts sequentially from top to bottom
- Read comments and text explanations carefully
- Modify code to deepen understanding
- Compare outputs with mathematical expectations

## Dependencies

- **MATLAB R2025a or newer** (for plain text Live Script support)
- **Symbolic Math Toolbox** (used in some scripts for symbolic computation)
- No additional toolboxes required for basic functionality

## Tips for Success

1. **Understand Before Coding:** Read the mathematical background before examining code
2. **Experiment:** Modify parameters to see how algorithms respond
3. **Visualize:** Pay attention to plots—they reveal algorithm behavior
4. **Compare:** Run multiple scripts and compare different approximation approaches
5. **Challenge Yourself:** Try to predict output before running code

## Syllabus Connection

This course emphasizes:
- Bridging theory and implementation
- Numerical stability and error analysis
- Algorithm development for engineering problems
- Scientific computing best practices

For complete course details, see `Syllabus_MECH.5200.pdf`

---

**Last Updated:** January 2026  
**Course:** MECH.5200 - Numerical Methods in Computational Mechanics
