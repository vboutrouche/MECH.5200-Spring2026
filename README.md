# MECH.5200 - Numerical Methods in Computational Mechanics

Welcome to the course materials for MECH.5200. This repository contains teaching scripts, code examples, and activities for learning numerical methods in computational mechanics.

## 📋 Assignments

### Project Proposal
**Due Date:** 13-Feb
**Document:** [Project.pdf](Misc/Project.pdf)

### Homework 1
**Due Date:** 4-Mar
**Document:** [Homework1.pdf](Misc/Homework1.pdf)

## 📅 Course Resources

- **[Class Schedule](Misc/ClassSchedule.pdf)** - Weekly topics, assignments, and important dates
- **[Course Syllabus](Misc/Syllabus_MECH.5200.pdf)** - Complete course information and policies

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
  - `Week2_Activity1_Solution.m` - Solutions for Activity 1
  - `Week2_Activity2.m` - Applications and comparisons of approximation methods

- **Video Supplements**
  - `Video3.pdf` - Supporting lecture material
  - `Video4.pdf` - Additional examples and derivations
  - `Video5.pdf` - Advanced topics and applications

### Week 3: [Content Title]
**Location:** `Week3/`

Additional course materials covering advanced topics in numerical methods.
- `Video6.pdf` - Lecture notes
- `Video7.pdf` - Additional materials
- `Video8.pdf` - Supplementary content


**Last Updated:** 4-Feb
**Course:** MECH.5200 - Numerical Methods in Computational Mechanics
