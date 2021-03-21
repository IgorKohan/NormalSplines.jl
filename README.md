# 1-D Normal Splines 

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://IgorKohan.github.io/NormalSplines.jl/stable)
[![Build Status](https://travis-ci.com/IgorKohan/NormalSplines.jl.svg?branch=master)](https://travis-ci.com/github/IgorKohan/NormalSplines.jl)
[![codecov.io](https://codecov.io/github/IgorKohan/NormalSplines.jl/coverage.svg?branch=master)](https://codecov.io/github/IgorKohan/NormalSplines.jl?branch=master)

This package implements the normal splines method for solving following interpolation problem:
![Problem definition](/images/problem-definition.gif)

## Instalation
Start Julia and run the following commands:
```
julia> using Pkg
julia> Pkg.add("NormalSplines")
```
## Example usage
To use the NormalSplines package, begin your code with
```
using NormalSplines
```
Construct a normal spline to some function values (knots ```x``` can be unevenly distributed):
```
using NormalSplines

x = [0.0, 1.0, 2.0]
u = [0.0, 1.0, 4.0]
interpolate(x, u, RK_W3())
```
here ```RK_W3()``` is the reproducing kernel of Sobolev space ```W_2^3[0,2]```.

Evaluate the spline, its first and second derivatives at some points:
```
p = [0.0, 0.5, 1.0, 1.5, 2.0]

σ = evaluate(p)       # result = [0.0, 0.26, 1.0, 2.24, 4.0]
σ = evaluate(1.5)     # result = 2.24

σ' = evaluate(p, 1)  # result = [0.08, 0.98, 1.96, 3.0, 4.04]

σ'' = evaluate(p, 2) # result = [1.71, 1.89, 2.02, 2.08, 2.08]
```

Construct a normal spline to some function and its first and second derivatives values:
```
using NormalSplines

x = [0.0, 1.0, 2.0] # Function knots
u = [0.0, 1.0, 4.0] # Function values 
s = [2.0]           # First derivative knot
v = [4.0]           # First derivative value
t = [0.0, 1.0]      # Second derivative knots
w = [2.0 ,2.0]      # Second derivative values
interpolate(x, u, s, v, t, w, RK_W3())
```
Evaluate the spline, its first and second derivatives at some points:
```
p = [0.0, 0.5, 1.0, 1.5, 2.0]

σ = evaluate(p)       # result = [0.0, 0.25, 1.0, 2.25, 4.0]

σ' = evaluate(p, 1) # result = [0.0, 1.0, 2.0, 3.0, 4.0]

σ'' = evaluate(p, 2) # result = [2.0, 2.0, 2.0, 2.0, 2.0]
```
Further examples are given in documentation.

## Documentation

For more information see [Documentation](https://igorkohan.github.io/NormalSplines.jl/stable).
