var documenterSearchIndex = {"docs":
[{"location":"Usage/#Example-Usage-1","page":"Example Usage","title":"Example Usage","text":"","category":"section"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"Let's interpolate function f(x)","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"f(x) =\n    begincases\n       0     1 le x lt 6 \n       1     6 le x le 10 \n       -x5 + 3     6 le x le 15 \n       0     15 le x le 20 \n    endcases","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"by values of the function in knots 1 2 3  20 (case A) and by values of the function and values of its first and second derivatives in the same knots (case B).","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"DocTestSetup = quote\n    using NormalSplines\nend","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"A)","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"    using NormalSplines\n\n    x = collect(1.0:1.0:20)       # function knots\n    u = x.*0.0                    # function values in knots\n    for i in 6:10\n        u[i] = 1.0\n    end\n    for i in 11:14\n        u[i] = -0.2 * i + 3.0\n    end\n\n    # build a twice-differentiable spline \n    # by values of function in knots\n    interpolate(x, u, RK_W3())\n\n    p = collect(1.0:0.2:20)        # evaluation points\n    σ = evaluate(p)\n    σ = nothing  # hide","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"(Image: Example 1A)","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"Evaluate the spline at some points:","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"    σ = evaluate([3.1, 8.1, 18.1])","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"B)","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"    using NormalSplines\n\n    x = collect(1.0:1.0:20)       # function knots\n    u = x.*0.0                    # function values in knots\n    for i in 6:10\n        u[i] = 1.0\n    end\n    for i in 11:14\n        u[i] = -0.2 * i + 3.0\n    end\n\n    s = x                         # function first derivative knots\n    v = x.*0.0                    # function first derivative values\n    for i in 11:14\n        v[i] = -0.2\n    end\n    t = x                         # function second derivative knots\n    w = x.*0.0                    # function second derivative values\n\n    # build a twice-differentiable spline by values of function,\n    # and values of its first and second derivatives in knots\n    interpolate(x, u, s, v, t, w, RK_W3())\n\n    p = collect(1.0:0.2:20)      # evaluation points\n    σ = evaluate(p)\n    σ = nothing  # hide","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"(Image: Example 1B)","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"Evaluate the spline at some points:","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"    σ = evaluate([3.1, 8.1, 18.1])","category":"page"},{"location":"Usage/#Q-and-A-1","page":"Example Usage","title":"Q & A","text":"","category":"section"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"Q1. Question: The call","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"interpolate(x, u, RK_H3())","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"cause the following error: PosDefException: matrix is not positive definite; Cholesky factorization failed. What is a reason of the error and how to resolve it?","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"A1. Answer: Creating a Bessel Potential kernel object with omitted parameter ε means that this paramter will be estimated during interpolating procedure execution. It might happen that estimated value of the ε is too small and corresponding  Gram matrix of linear system of equations which defines the normal spline coefficients is very ill-conditioned and it lost its positive definiteness because of floating-point rounding errors.  ","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"There are two ways to fix it.","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"We can get the estimated value of ε by calling function get_epsilon():","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"ε = get_epsilon()","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"then we could try to call the interpolate function with greater ε value of the reproducing kernel parameter:","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"interpolate(x, u, RK_H3(10.0*ε))","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"We may change the precision of floating point calculations. Namely it is possible to use Julia standard BigFloat numbers or Double64 - extended precision float type from the package DoubleFloats:","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"using DoubleFloats\n\nx = Double64.(x)\nu = Double64.(u)\ninterpolate(x, u, RK_H3())","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"This answer also applies to types RK_H1() and RK_H2().","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"Q2. Question: The call","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"interpolate(x, u, RK_W3())","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"cause the following error: PosDefException: matrix is not positive definite; Cholesky factorization failed. How to resolve it?","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"A2. Answer: The reason of that exception is the Gram matrix of linear system of equations which defines the normal spline coefficients is very ill-conditioned and it lost its positive definiteness because of floating-point rounding errors.  ","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"The only way to fix it - is using floating-point arithmetic with extended precision. It can be provided by Julia standard BigFloat type or Double64 type from the package DoubleFloats:","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"using DoubleFloats\n\nx = Double64.(x)\nu = Double64.(u)\ninterpolate(x, u, RK_W3())","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"This answer also applies to types RK_W1() and RK_W2().","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"Q3. Question: The following calls","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"interpolate(x, u, RK_H3())\nσ = evaluate(p)","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"produce the output which is not quite satisfactoty. Is it possible to improve the quality of interpolation?","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"A3. Answer: Creating a Bessel Potential kernel object with omitted parameter ε means that this paramter will be estimated during interpolating procedure execution. It might happen that estimated value of the ε is too large and it is possible to use a smaller ε value which would result in better quality of interpolation. We can get the estimated value of ε by calling function get_epsilon():","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"ε = get_epsilon()","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"and get an estimation of the problem's Gram matrix condition number by calling get_cond() function:","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"cond = get_cond()","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"In a case when estimated condition number is not very large, i.e. less than 10^12 using standard Float64 floating-point arithmetic, we may attempt to build a better interpolation spline by calls:","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"interpolate(x, u, RK_H3(ε/5))\nσ = evaluate(p)","category":"page"},{"location":"Usage/#","page":"Example Usage","title":"Example Usage","text":"Another option is using a smaller value of the ε and perform calculations using extended precision floating-point arithmetic.","category":"page"},{"location":"#NormalSplines.jl-1","page":"Home","title":"NormalSplines.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"1-D Normal Splines implementation in Julia","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\n      \"index.md\",\n      \"Usage.md\",\n      \"Public-API.md\",\n      \"Interpolating-Normal-Splines.md\",\n      \"Reproducing-Kernel-of-Bessel-Potential-space.md\"\n]\nDepth = 3","category":"page"},{"location":"#","page":"Home","title":"Home","text":"NormalSplines.jl implements the normal splines method for solving following interpolation problem:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Problem: Given points x_1 lt  x_2 lt dots lt x_n_1, s_1 lt s_2 lt dots lt s_n_2 and  t_1 lt t_2 lt dots lt t_n_3 find a function f such that","category":"page"},{"location":"#","page":"Home","title":"Home","text":"tag1\nbeginaligned\n f(x_i) =  u_i   quad  i = 1 2 dots n_1  \n  \n f(s_j) =  v_j   quad  j = 1 2 dots n_2  \n  \n f(t_k) =  w_k   quad  k = 1 2 dots n_3  \n\n n_1 gt 0       n_2 ge 0      n_3 ge 0  \nendaligned","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Knots x_i, s_j and t_k may coincide. We assume that function f is an element of an appropriate reproducing kernel Hilbert space H. ","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The normal splines method consists in finding a solution of system (1) having minimal norm in Hilbert space H, thus the interpolation normal spline sigma is defined as follows:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"tag2\n   sigma = rm argmin   f ^2_H  (1) forall f in H   ","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Normal splines method is based on the following functional analysis results:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Embedding theorems (Sobolev embedding theorem and Bessel potential space embedding theorem)\nThe Riesz representation theorem for Hilbert spaces\nReproducing kernel properties. ","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Using these results it is possible to reduce task (2) to solving a system of linear equations with symmetric positive definite Gram matrix.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Normal splines are constructed in Sobolev space W^l_2 a b with norm:","category":"page"},{"location":"#","page":"Home","title":"Home","text":" f _W_2^l = left ( sum_i=0^l-1 (f^(i)(a))^2 + int_a^b (f^(l)(s))^2 ds right )^12   quad l = 3 2 textor  1","category":"page"},{"location":"#","page":"Home","title":"Home","text":"and in Bessel potential space H_varepsilon^l(R) defined as","category":"page"},{"location":"#","page":"Home","title":"Home","text":"   H^l_varepsilon (R) = left f  f in S \n  ( varepsilon ^2 +  s ^2 )^l2mathcal F f in L_2 (R) right    varepsilon gt 0  quad l = 3 2 textor  1","category":"page"},{"location":"#","page":"Home","title":"Home","text":"where S  (R) is space of L. Schwartz tempered distributions and mathcal F f is a Fourier transform of the f. Space H^l_varepsilon (R) is Hilbert space with norm ","category":"page"},{"location":"#","page":"Home","title":"Home","text":" f _ H^l_varepsilon =   (  varepsilon ^2 +  s ^2 )^l2 mathcal F varphi  _L_2  ","category":"page"},{"location":"#","page":"Home","title":"Home","text":"If n_3  0 then value of l is 3, in a case when n_3 = 0 and n_2  0 the value of l must be 2 or 3.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Detailed explanation is given in Interpolating Normal Splines","category":"page"},{"location":"Public-API/#Public-API-1","page":"Public API","title":"Public API","text":"","category":"section"},{"location":"Public-API/#","page":"Public API","title":"Public API","text":"Order   = [:function, :type]","category":"page"},{"location":"Public-API/#Functions-1","page":"Public API","title":"Functions","text":"","category":"section"},{"location":"Public-API/#","page":"Public API","title":"Public API","text":"interpolate\nevaluate\nget_cond\nget_epsilon","category":"page"},{"location":"Public-API/#NormalSplines.interpolate","page":"Public API","title":"NormalSplines.interpolate","text":"interpolate(knots::Vector{T}, values::Vector{T}, kernel :: RK, factorize = true) where {T <: AbstractFloat, RK <: ReproducingKernel_1}\n\nCreate an interpolating normal spline by values of function in knots.\n\nArguments\n\nknots::Vector{T}: locations of given function values.\nvalues::Vector{T}: function values in knots.\nkernel::RK: reproducing kernel of Hilbert space the normal spline is constructed in.               It must be a struct object of the following type:                 RK_W1 or RK_H1 if the spline is constructing as a continuous function,                 RK_W2 or RK_H2 if the spline is constructing as a differentiable function,                 RK_W3 or RK_H3 if the spline is constructing as a twice differentiable function.\nfactorize::Bool: defines if it is necessary to compute the Gram matrix factorization.                 Must be set to true if the spline is constructing first time.                 Can be set to false if the spline is constructing with the same knots                 as it was done first time but with different function values in knots.\n\nReturns : nothing.\n\n\n\n\n\ninterpolate(knots::Vector{T}, values::Vector{T},              d_knots::Vector{T}, d_values::Vector{T},              kernel :: RK, derivative = 1, factorize = true) where {T <: AbstractFloat, RK <: ReproducingKernel_2}\n\nCreate an interpolating normal spline by values of function and it's first or second derivatives in knots.\n\nArguments\n\nknots::Vector{T}: locations of given function values.\nvalues::Vector{T}: function values in knots.\nd_knots::Vector{T}: locations of given function first derivative values (if derivative=1)                       or second derivative values (if derivative=2).\nd_values::Vector{T}: values of function first derivative (if derivative=1) or                        values of function second derivative (if derivative=2) in d_knots.\nkernel::RK: reproducing kernel of Hilbert space the normal spline is constructed in.               In a case when derivative=1, kernel must be a struct object of the following type:                    RK_W2 or RK_H2 if the spline is constructing as a differentiable function,                    RK_W3 or RK_H3 if the spline is constructing as a twice differentiable function.               In a case when derivative=2, kernel must be a struct object of the following type:                    RK_W3 or RK_H3 - the spline is constructing as a twice differentiable function,\nderivative::Int: 1 - d_knots and d_values provide data of function first derivative,                    2 - d_knots and d_values provide data of function second derivative.\nfactorize::Bool: defines if it is necessary to compute the Gram matrix factorization.                    Must be set to true if the spline is constructing first time.                    Can be set to false if the spline is constructing with the same knots                    as it was done first time but with different function values in knots.\n\nReturns : nothing.\n\n\n\n\n\ninterpolate(knots::Vector{T}, values::Vector{T},              d1_knots::Vector{T}, d1_values::Vector{T},              d2_knots::Vector{T}, d2_values::Vector{T},              kernel :: RK, factorize = true) where {T <: AbstractFloat, RK <: ReproducingKernel_3}\n\nCreate an interpolating normal spline by values of function and it's first and second derivatives in knots.\n\nArguments\n\nknots::Vector{T}: locations of given function values.\nvalues::Vector{T}: function values in knots.\nd1_knots::Vector{T}: locations of given function derivative values.\nd1_values::Vector{T}: values of function derivative in d1_knots.\nd2_knots::Vector{T}: locations of given function second derivative values.\nd2_values::Vector{T}: values of function second derivative in d2_knots.\nkernel::RK: reproducing kernel of Hilbert space the normal spline is constructed in.               It must be a struct object of the following type:                 RK_W3 or RK_H3 - the spline is constructing as a twice differentiable function.\nfactorize::Bool: defines if it is necessary to compute the Gram matrix factorization.                 Must be set to true if the spline is constructing first time.                 Can be set to false if the spline is constructing with the same knots                 as it was done first time but with different function values in knots.\n\nReturns : nothing.\n\n\n\n\n\n","category":"function"},{"location":"Public-API/#NormalSplines.evaluate","page":"Public API","title":"NormalSplines.evaluate","text":"evaluate(points::Vector{T}, derivative::Int = 0) where T <: AbstractFloat\n\nEvaluate normal spline and it's derivatives.\n\nArguments\n\npoints::Vector{T}: locations at which spline or its derivative are evaluating.\nderivative::Int: 0 - spline evaluation.                    1 - first derivative of spline evaluation.                    2 - second derivative of spline evaluation.\n\nReturns: Vector{T} of spline values or its derivatives at the locations           defined in points.\n\n\n\n\n\nevaluate(point::T, derivative::Int = 0) where T <: AbstractFloat\n\nEvaluate normal spline and it's derivatives.\n\nArguments\n\npoint::T: location at which spline or its derivative are evaluating.\nderivative::Int: 0 - spline evaluation.                    1 - first derivative of spline evaluation.                    2 - second derivative of spline evaluation.\n\nReturns: spline value or its derivative at the point location.\n\n\n\n\n\n","category":"function"},{"location":"Public-API/#NormalSplines.get_cond","page":"Public API","title":"NormalSplines.get_cond","text":"get_cond()\n\nEstimate condition number of Gram matrix.\n\n\n\n\n\n","category":"function"},{"location":"Public-API/#NormalSplines.get_epsilon","page":"Public API","title":"NormalSplines.get_epsilon","text":"get_epsilon()\n\nGet value of parameter ε of struct object RK_H1, RK_H2 or RK_H3.\n\n\n\n\n\n","category":"function"},{"location":"Public-API/#Types-1","page":"Public API","title":"Types","text":"","category":"section"},{"location":"Public-API/#Sobolev-space-Reproducing-Kernels-1","page":"Public API","title":"Sobolev space Reproducing Kernels","text":"","category":"section"},{"location":"Public-API/#","page":"Public API","title":"Public API","text":"RK_W1\nRK_W2\nRK_W3","category":"page"},{"location":"Public-API/#NormalSplines.RK_W1","page":"Public API","title":"NormalSplines.RK_W1","text":"struct RK_W1 <: ReproducingKernel_1\n\nDefine a type of reproducing kernel of Sobolev space W^2_1 01:\n\nV(eta xi) =\n  begincases\n       1 + eta      0 le eta le xi le 1 \n       1 + xi      0 le xi le eta le 1\n    endcases\n\n\n\n\n\n","category":"type"},{"location":"Public-API/#NormalSplines.RK_W2","page":"Public API","title":"NormalSplines.RK_W2","text":"struct RK_W2 <: ReproducingKernel_2\n\nDefine a type of reproducing kernel of Sobolev space W^2_2 01:\n\nV(eta xi) =\n  begincases\n       1 + (eta + eta^2  2)xi - eta^3  6       0 le eta le xi le 1  \n       1 + (xi + xi^2  2)eta - xi^3  6      0 le xi le eta le 1\n    endcases\n\n\n\n\n\n","category":"type"},{"location":"Public-API/#NormalSplines.RK_W3","page":"Public API","title":"NormalSplines.RK_W3","text":"struct RK_W3 <: ReproducingKernel_3\n\nDefine a type of reproducing kernel of Sobolev space W^2_3 01:\n\nV(eta xi) =\n  begincases\n       sum_i=0^2 fracxi^ii left ( fraceta^ii + (-1)^i frac eta^5 - i(5 - i)  right )       0 le eta le xi le 1\n        \n       sum_i=0^2 fraceta^ii left ( fracxi^ii + (-1)^i frac xi^5 - i(5 - i)  right )      0 le xi le eta le 1   endcases\n\n\n\n\n\n","category":"type"},{"location":"Public-API/#Bessel-Potential-space-Reproducing-Kernels-1","page":"Public API","title":"Bessel Potential space Reproducing Kernels","text":"","category":"section"},{"location":"Public-API/#","page":"Public API","title":"Public API","text":"RK_H1\nRK_H2\nRK_H3","category":"page"},{"location":"Public-API/#NormalSplines.RK_H1","page":"Public API","title":"NormalSplines.RK_H1","text":"struct RK_H1{T <: AbstractFloat} <: ReproducingKernel_1\n\nDefine a type of reproducing kernel of Bessel Potential space H^1_ε (R):\n\nV(eta  xi varepsilon) = exp (-varepsilon xi - eta)  \n\nArguments\n\nno argument: means that parameter ε will be estimated in course of              interpolation procedure execution.\nε::T: 'scaling parameter' from the Bessel Potential space definition,          it must be greater than zero.\n\n\n\n\n\n","category":"type"},{"location":"Public-API/#NormalSplines.RK_H2","page":"Public API","title":"NormalSplines.RK_H2","text":"struct RK_H2{T <: AbstractFloat} <: ReproducingKernel_2\n\nDefine a type of reproducing kernel of Bessel Potential space H^2_ε (R):\n\nV(eta  xi varepsilon) = exp (-varepsilon xi - eta)\n             (1 + varepsilon xi  - eta)  \n\nArguments\n\nno argument: means that parameter ε will be estimated in course of              interpolation procedure execution.\nε::T: 'scaling parameter' from the Bessel Potential space definition,          it must be greater than zero.\n\n\n\n\n\n","category":"type"},{"location":"Public-API/#NormalSplines.RK_H3","page":"Public API","title":"NormalSplines.RK_H3","text":"struct RK_H3{T <: AbstractFloat} <: ReproducingKernel_3\n\nDefine a type of reproducing kernel of Bessel Potential space H^3_ε (R):\n\nV(eta  xi varepsilon) = exp (-varepsilon xi - eta)\n             (3 + 3varepsilon xi  - eta + varepsilon ^2 xi - eta ^2 )  \n\nArguments\n\nno argument: means that parameter ε will be estimated in course of              interpolation procedure execution.\nε::T: 'scaling parameter' from the Bessel Potential space definition,          it must be greater than zero.\n\n\n\n\n\n","category":"type"}]
}
