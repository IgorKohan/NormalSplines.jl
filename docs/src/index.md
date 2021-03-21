
# NormalSplines.jl

*1-D Normal Splines implementation in Julia*

```@contents
Pages = [
      "index.md",
      "Usage.md",
      "Public-API.md",
      "Interpolating-Normal-Splines.md",
      "Reproducing-Kernel-of-Bessel-Potential-space.md"
]
Depth = 3
```

`NormalSplines.jl` implements the normal splines method for solving following interpolation problem:

*Problem:* Given points ``x_1 \lt  x_2 \lt \dots \lt x_{n_1}``, ``s_1 \lt s_2 \lt \dots \lt s_{n_2}`` and  ``t_1 \lt t_2 \lt \dots \lt t_{n_3}`` find a function ``f`` such that

```math
\tag{1}
\begin{aligned}
& f(x_i) =  u_i \, , \quad  i = 1, 2, \dots, n_1 \, ,
\\  
& f'(s_j) =  v_j \, , \quad  j = 1, 2, \dots, n_2 \, ,
\\  
& f''(t_k) =  w_k \, , \quad  k = 1, 2, \dots, n_3 \, ,
\\
& n_1 \gt 0 \, ,  \ \  n_2 \ge 0 \, , \ \  n_3 \ge 0 \, .
\end{aligned}
```
Knots ``\{x_i\}``, ``\{s_j\}`` and ``\{t_k\}`` may coincide. We assume that function ``f`` is an element of an appropriate reproducing kernel Hilbert space ``H``. 

The normal splines method consists in finding a solution of system (1) having minimal norm in Hilbert space ``H``, thus the interpolation normal spline ``\sigma`` is defined as follows:

```math
\tag{2}
   \sigma = {\rm arg\,min}\{  \| f \|^2_H : (1), \forall f \in H \} \, .
```

Normal splines method is based on the following functional analysis results:

* Embedding theorems (Sobolev embedding theorem and Bessel potential space embedding theorem)
* The Riesz representation theorem for Hilbert spaces
* Reproducing kernel properties. 

Using these results it is possible to reduce task (2) to solving a system of linear equations with symmetric positive definite Gram matrix.

Normal splines are constructed in Sobolev space ``W^l_2 [a, b]`` with norm:

```math
\| f \|_{W_2^l} = \left ( \sum_{i=0}^{l-1} (f^{(i)}(a))^2 + \int_a^b (f^{(l)}(s))^2 ds \right )^{1/2} \, , \quad l = 3, 2, \text{or} \ 1,
```
and in Bessel potential space ``H_\varepsilon^l(R)`` defined as

```math
   H^l_\varepsilon (R) = \left\{ f | f \in S' ,
  ( \varepsilon ^2 + | s |^2 )^{l/2}{\mathcal F} [f] \in L_2 (R) \right\} \, ,  \varepsilon \gt 0 , \quad l = 3, 2, \text{or} \ 1,
```
where ``S'  (R)`` is space of L. Schwartz tempered distributions and ``\mathcal F [f]`` is a Fourier transform of the ``f``. Space ``H^l_\varepsilon (R)`` is Hilbert space with norm 

```math
\| f \|_ {H^l_\varepsilon} =  \| (  \varepsilon ^2 + | s |^2 )^{l/2} {\mathcal F} [\varphi ] \|_{L_2} \ .
```
If ``n_3 > 0`` then value of ``l`` is ``3``, in a case when ``n_3 = 0`` and ``n_2 > 0`` the value of ``l`` must be ``2`` or ``3``.

Detailed explanation is given in [Interpolating Normal Splines](https://igorkohan.github.io/NormalSplines.jl/stable/Interpolating-Normal-Splines/)
