_kernel = nothing
_ε = nothing
_knots = nothing
_d1_knots = nothing
_d2_knots = nothing
_mat = nothing
_chol = nothing
_mu = nothing
_mu_1 = nothing
_mu_2 = nothing

_n = 0
_n_1 = 0
_n_2 = 0

_isOK = false
_isPositiveDefinite = true
_cond = 0.0

_a = 0.0
_b = 0.0

function _cleanup()
    global _a = 0.0
    global _b = 0.0
    global _n = 0
    global _n_1 = 0
    global _n_2 = 0
    global _kernel = nothing
    global _ε = nothing
    global _knots = nothing
    global _d1_knots = nothing
    global _d2_knots = nothing
    global _chol = nothing
    global _isOK = false
    global _isPositiveDefinite = true
    global _cond = 0.0
    global _mu = nothing
    global _mu_1 = nothing
    global _mu_2 = nothing
end

function _estimate_ε(knots::Vector{T}) where T <: AbstractFloat
    s = T(0.0)
    n = length(knots)
    for i = 1:n
        for j = i:n
            s += abs(knots[i] - knots[j])
        end
    end
    if s > T(0.0)
        s *= T(1.0) / (T(n)* sqrt(T(n)))
    else
        s = T(1.0)
    end
    global _ε = s
    return s
end

function _get_epsilon()
    return _ε
end

function _prepare_spline(kernel::RK,
                         knots::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
    _cleanup()

    a = minimum(knots)
    b = maximum(knots)
    d = b - a
    n = length(knots)
    kns = Vector{T}(undef, n)
    for k = 1:n
        kns[k] = (knots[k] - a)/d
    end
    if isa(kernel, RK_H1) || isa(kernel, RK_H2) || isa(kernel, RK_H3)
        global _ε = kernel.ε
        if T(kernel.ε) == T(0.0)
            ε = _estimate_ε(kns)
            if isa(kernel, RK_H1)
                kernel = RK_H1(ε)
            end
            if isa(kernel, RK_H2)
                kernel = RK_H2(ε)
            end
            if isa(kernel, RK_H3)
                kernel = RK_H3(ε)
            end
        end
    end
    mat = _gram(kernel, kns)
    _cholesky_factorization(mat)
    global _kernel = kernel
    global _mat = mat
    global _knots = kns
    global _n = n
    global _a = a
    global _b = b
    return
end

function _prepare_spline_d1(kernel::RK,
                            knots::Vector{T},
                            d1_knots::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_2}
    _cleanup()

    if isa(kernel, RK_H1) || isa(kernel, RK_W1)
       Throw(ArgumentError("kernel: Incorrect parameter type (must be one of RK_H2, RK_H3, RK_W2, RK_W3)."))
    end
    a = minimum(knots)
    b = maximum(knots)
    a_1 = minimum(d1_knots)
    b_1 = maximum(d1_knots)
    a = min(a, a_1)
    b = max(b, b_1)
    d = b - a
    n = length(knots)
    kns = Vector{T}(undef, n)
    for k = 1:n
        kns[k] = (knots[k] - a)/d
    end
    n_1 = length(d1_knots)
    d1_kns = Vector{T}(undef, n_1)
    for k = 1:n_1
        d1_kns[k] = (d1_knots[k] - a)/d
    end
    if isa(kernel, RK_H2) || isa(kernel, RK_H3)
        global _ε = kernel.ε
        if T(kernel.ε) == T(0.0)
            ε = _estimate_ε(kns)
            if isa(kernel, RK_H2)
                kernel = RK_H2(ε)
            end
            if isa(kernel, RK_H3)
                kernel = RK_H3(ε)
            end
        end
    end
    mat = _gram_d1(kernel, kns, d1_kns)
    _cholesky_factorization(mat)
    global _kernel = kernel
    global _mat = mat
    global _knots = kns
    global _d1_knots = d1_kns
    global _n = n
    global _n_1 = n_1
    global _a = a
    global _b = b
    return
end

function _prepare_spline_d2(kernel::RK,
                            knots::Vector{T},
                            d2_knots::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
    _cleanup()
    if isa(kernel, RK_H1) || isa(kernel, RK_W1) || isa(kernel, RK_H2) || isa(kernel, RK_W2)
       Throw(ArgumentError("kernel: Incorrect parameter type (must be RK_H3 or RK_W3)."))
    end
    a = minimum(knots)
    b = maximum(knots)
    a_2 = minimum(d2_knots)
    b_2 = maximum(d2_knots)
    a = min(a, a_2)
    b = max(b, b_2)
    d = b - a
    n = length(knots)
    kns = Vector{T}(undef, n)
    for k = 1:n
        kns[k] = (knots[k] - a)/d
    end
    n_2 = length(d2_knots)
    d2_kns = Vector{T}(undef, n_2)
    for k = 1:n_2
        d2_kns[k] = (d2_knots[k] - a)/d
    end
    if isa(kernel, RK_H3)
        global _ε = kernel.ε
        if T(kernel.ε) == T(0.0)
            ε = _estimate_ε(kns)
            kernel = RK_H3(ε)
        end
    end
    mat = _gram_d2(kernel, kns, d2_kns)
    _cholesky_factorization(mat)
    global _kernel = kernel
    global _mat = mat
    global _knots = kns
    global _d2_knots = d2_kns
    global _n = n
    global _n_2 = n_2
    global _a = a
    global _b = b
    return
end

function _prepare_spline_d12(kernel::RK,
                             knots::Vector{T},
                             d1_knots::Vector{T},
                             d2_knots::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
    _cleanup()
    if isa(kernel, RK_H1) || isa(kernel, RK_W1) || isa(kernel, RK_H2) || isa(kernel, RK_W2)
       Throw(ArgumentError("kernel: Incorrect parameter type (must be RK_H3 or RK_W3)."))
    end
    a = minimum(knots)
    b = maximum(knots)
    a_1 = minimum(d1_knots)
    b_1 = maximum(d1_knots)
    a_2 = minimum(d2_knots)
    b_2 = maximum(d2_knots)
    a = min(a, a_1, a_2)
    b = max(b, b_1, b_2)
    d = b - a
    n = length(knots)
    kns = Vector{T}(undef, n)
    for k = 1:n
        kns[k] = (knots[k] - a)/d
    end
    n_1 = length(d1_knots)
    d1_kns = Vector{T}(undef, n_1)
    for k = 1:n_1
        d1_kns[k] = (d1_knots[k] - a)/d
    end
    n_2 = length(d2_knots)
    d2_kns = Vector{T}(undef, n_2)
    for k = 1:n_2
        d2_kns[k] = (d2_knots[k] - a)/d
    end
    if isa(kernel, RK_H3)
        global _ε = kernel.ε
        if T(kernel.ε) == T(0.0)
            ε = _estimate_ε(kns)
            kernel = RK_H3(ε)
        end
    end
    mat = _gram12(kernel, kns, d1_kns, d2_kns)
    _cholesky_factorization(mat)
    global _kernel = kernel
    global _mat = mat
    global _knots = kns
    global _d1_knots = d1_kns
    global _d2_knots = d2_kns
    global _n = n
    global _n_1 = n_1
    global _n_2 = n_2
    global _a = a
    global _b = b
    return
end

function _construct_spline(values::Vector{T}) where T <: AbstractFloat
    if(length(values) != _n)
        throw(DimensionMismatch("Number of data values does not correspond to the number of knots."))
    end
    global _mu = nothing
    if !_isOK
        if !_isPositiveDefinite
            throw(LinearAlgebra.PosDefException(0))
            return
        end
        error("Gram matrix was not factorized.")
        return
    end
    chol = _chol
    mu = Vector{T}(undef, size(chol, 1))
    ldiv!(mu, chol, values)
    global _mu = mu
    return
end

function _construct_spline_d1(values::Vector{T},
                              d1_values::Vector{T}) where T <: AbstractFloat
    if(length(values) != _n)
        throw(DimensionMismatch("Number of data values does not correspond to the number of knots."))
    end
    if(length(d1_values) != _n_1)
        throw(DimensionMismatch("Number of data d1-values does not correspond to the number of d1-knots."))
    end
    global _mu = nothing
    global _mu_1 = nothing
    global _mu_2 = nothing
    if !_isOK
        if !_isPositiveDefinite
            throw(LinearAlgebra.PosDefException(0))
            return
        end
        error("Gram matrix was not factorized.")
        return
    end

    n_1 = _n_1
    d = _b - _a
    chol = _chol
    muf = Vector{T}(undef, size(chol, 1))
    b = Vector{T}(values)
    d1_vals = Vector{T}(d1_values)
    d1_vals .*= d
    append!(b, d1_vals)
    ldiv!(muf, chol, b)

    global _mu = muf[1:_n]
    global _mu_1 = muf[(_n + 1):end]
    return
end

function _construct_spline_d2(values::Vector{T},
                              d2_values::Vector{T}) where T <: AbstractFloat
    if(length(values) != _n)
        throw(DimensionMismatch("Number of data values does not correspond to the number of knots."))
    end
    if(length(d2_values) != _n_2)
        throw(DimensionMismatch("Number of data d2-values does not correspond to the number of d2-knots."))
    end
    global _mu = nothing
    global _mu_1 = nothing
    global _mu_2 = nothing
    if !_isOK
        if !_isPositiveDefinite
            throw(LinearAlgebra.PosDefException(0))
            return
        end
        error("Gram matrix was not factorized.")
        return
    end

    n_2 = _n_2
    d = _b - _a
    d2 = d^2

    chol = _chol
    muf = Vector{T}(undef, size(chol, 1))
    b = Vector{T}(values)
    d2_vals = Vector{T}(d2_values)
    d2_vals .*= d2
    append!(b, d2_vals)
    ldiv!(muf, chol, b)

    global _mu = muf[1:_n]
    global _mu_2 = muf[(_n + 1):end]
    return
end

function _construct_spline_d12(values::Vector{T},
                               d1_values::Vector{T},
                               d2_values::Vector{T}) where T <: AbstractFloat
    if(length(values) != _n)
        throw(DimensionMismatch("Number of data values does not correspond to the number of knots."))
    end
    if(length(d1_values) != _n_1)
        throw(DimensionMismatch("Number of data d1-values does not correspond to the number of d1-knots."))
    end
    if(length(d2_values) != _n_2)
        throw(DimensionMismatch("Number of data d2-values does not correspond to the number of d2-knots."))
    end
    global _mu = nothing
    global _mu_1 = nothing
    global _mu_2 = nothing
    if !_isOK
        if !_isPositiveDefinite
            throw(LinearAlgebra.PosDefException(0))
            return
        end
        error("Gram matrix was not factorized.")
        return
    end

    n_1 = _n_1
    d = _b - _a

    n_2 = _n_2
    d2 = d^2

    chol = _chol
    muf = Vector{T}(undef, size(chol, 1))
    b = Vector{T}(values)
    d1_vals = Vector{T}(d1_values)
    d1_vals .*= d
    d2_vals = Vector{T}(d2_values)
    d2_vals .*= d2
    append!(b, d1_vals)
    append!(b, d2_vals)
    ldiv!(muf, chol, b)
    global _mu = muf[1:_n]
    global _mu_1 = muf[(_n + 1):(_n + _n_1)]
    global _mu_2 = muf[(_n + _n_1 + 1):end]
    return
end

function _evaluate(points::Vector{T}) where T <: AbstractFloat
    if _mu == nothing
        error("Spline coefficients were not calculated.")
    end

    kernel = deepcopy(_kernel)
    mu = Vector{T}(_mu)
    n = length(mu)
    n_1 = _mu_1 != nothing ? length(_mu_1) : 0
    n_2 = _mu_2 != nothing ? length(_mu_2) : 0

    h_values = Vector{T}(undef, n)
    knots = Vector{T}(_knots)
    m = length(points)
    spline_values = Vector{T}(undef, m)
    pnts = copy(points)
    a = _a
    b = _b
    d = b - a
    for k = 1:m
        pnts[k] = (pnts[k] - a)/d
    end
    @inbounds for p = 1:m
        for i = 1:n
            h_values[i] = _rk(kernel, pnts[p], knots[i])
        end
        spline_values[p] = mu ⋅ h_values
    end
    if n_1 > 0
        d1_knots = Vector{T}(_d1_knots)
        mu_1 = Vector{T}(_mu_1)
        h1_values = Vector{T}(undef, n_1)
        @inbounds for p = 1:m
            for i = 1:n_1
                h1_values[i] = _∂rk_∂ξ(kernel, pnts[p], d1_knots[i])
            end
            spline_values[p] += mu_1 ⋅ h1_values
        end
    end
    if n_2 > 0
        d2_knots = Vector{T}(_d2_knots)
        mu_2 = Vector{T}(_mu_2)
        h2_values = Vector{T}(undef, n_2)
        @inbounds for p = 1:m
            for i = 1:n_2
                h2_values[i] = _∂²rk_∂ξ²(kernel, pnts[p], d2_knots[i])
            end
            spline_values[p] += mu_2 ⋅ h2_values
        end
    end
    return spline_values
end

function _evaluate_d1(points::Vector{T}) where T <: AbstractFloat
    if _mu == nothing
        error("Spline coefficients were not calculated.")
    end
    if isa(_kernel, RK_W1) || isa(_kernel, RK_H1)
        error("Spline was not built as a differentiable one.")
    end

    kernel = deepcopy(_kernel)
    mu = Vector{T}(_mu)
    n = length(mu)
    n_1 = _mu_1 != nothing ? length(_mu_1) : 0
    n_2 = _mu_2 != nothing ? length(_mu_2) : 0

    m = length(points)
    spline_values = Vector{T}(undef, m)
    h_values = Vector{T}(undef, n)
    knots = Vector{T}(_knots)
    pnts = copy(points)
    a = _a
    b = _b
    d = b - a
    for k = 1:m
        pnts[k] = (pnts[k] - a)/d
    end
    @inbounds for p = 1:m
        for i = 1:n
            h_values[i] = _∂rk_∂η(kernel, pnts[p], knots[i])
        end
        spline_values[p] = mu ⋅ h_values
    end
    if n_1 > 0
        d1_knots = Vector{T}(_d1_knots)
        mu_1 = Vector{T}(_mu_1)
        h1_values = Vector{T}(undef, n_1)
        @inbounds for p = 1:m
            for i = 1:n_1
                h1_values[i] = _∂²rk_∂η∂ξ(kernel, pnts[p], d1_knots[i])
            end
            spline_values[p] += mu_1 ⋅ h1_values
        end
    end
    if n_2 > 0
        d2_knots = Vector{T}(_d2_knots)
        mu_2 = Vector{T}(_mu_2)
        h2_values = Vector{T}(undef, n_2)
        @inbounds for p = 1:m
            for i = 1:n_2
                h2_values[i] = _∂³rk_∂η∂ξ²(kernel, pnts[p], d2_knots[i])
            end
            spline_values[p] += mu_2 ⋅ h2_values
        end
    end

    spline_values ./= d
    return spline_values
end

function _evaluate_d2(points::Vector{T}) where T <: AbstractFloat
    if _mu == nothing
        error("Spline coefficients were not calculated.")
    end
    if !isa(_kernel, RK_W3) && !isa(_kernel, RK_H3)
        error("Spline was not built as a twice differentiable one.")
    end

    kernel = deepcopy(_kernel)
    mu = Vector{T}(_mu)
    n = length(mu)
    n_1 = _mu_1 != nothing ? length(_mu_1) : 0
    n_2 = _mu_2 != nothing ? length(_mu_2) : 0

    m = length(points)
    spline_values = Vector{T}(undef, m)
    h_values = Vector{T}(undef, n)
    knots = Vector{T}(_knots)
    pnts = copy(points)
    a = _a
    b = _b
    d = b - a
    for k = 1:m
        pnts[k] = (pnts[k] - a)/d
    end
    @inbounds for p = 1:m
        for i = 1:n
            h_values[i] = _∂²rk_∂η²(kernel, pnts[p], knots[i])
        end
        spline_values[p] = mu ⋅ h_values
    end
    if n_1 > 0
        d1_knots = Vector{T}(_d1_knots)
        mu_1 = Vector{T}(_mu_1)
        h1_values = Vector{T}(undef, n_1)
        @inbounds for p = 1:m
            for i = 1:n_1
                h1_values[i] = _∂³rk_∂η²∂ξ(kernel, pnts[p], d1_knots[i])
            end
            spline_values[p] += mu_1 ⋅ h1_values
        end
    end
    if n_2 > 0
        d2_knots = Vector{T}(_d2_knots)
        mu_2 = Vector{T}(_mu_2)
        h2_values = Vector{T}(undef, n_2)
        @inbounds for p = 1:m
            for i = 1:n_2
                h2_values[i] = _∂⁴rk_∂η²∂ξ²(kernel, pnts[p], d2_knots[i])
            end
            spline_values[p] += mu_2 ⋅ h2_values
        end
    end

    d2 = d^2
    spline_values ./= d2
    return spline_values
end

# ```
# Get estimation of the Gram matrix condition number
# Brás, C.P., Hager, W.W. & Júdice, J.J. An investigation of feasible descent algorithms for estimating the condition number of a matrix. TOP 20, 791–809 (2012).
# https://link.springer.com/article/10.1007/s11750-010-0161-9
# ```
function _get_cond(nit = 3)
    cond =  0.0
    if !_isOK
        error("_get_cond: Gram matrix was not factorized.")
    end
    mat = _mat
    T = eltype(mat)
    mat_norm = norm(mat, 1)
    chol = _chol
    n = size(mat, 1)
    x = Vector{T}(undef, n)
    x .= 1 / T(n)
    z = Vector{T}(undef, n)
    gamma = T(0)
    for it = 1:nit
        z = ldiv!(z, chol, x)
        gamma = T(0);
        for i = 1:n
            gamma += abs(z[i])
            z[i] = sign(z[i])
        end
        z = ldiv!(z, chol, copy(z))
        zx = z ⋅ x
        idx = 1
        for i = 1:n
            z[i] = abs(z[i])
            if z[i] > z[idx]
                idx = i
            end
        end
        if z[idx] <= zx
            break
        end
        x .= T(0)
        x[idx] = T(1);
    end
    cond = 10.0^floor(log10(mat_norm * gamma));
    global _cond = cond
    return cond
end

function _get_result_precision()
    T = eltype(_mat)
    cond = _get_cond()
    precision = -Int64(floor(log10(eps(T) * cond)))
    return precision
end

function _cholesky_factorization(mat::Matrix{T}) where T <: AbstractFloat
    try
        chol = cholesky(mat)
        global _isOK = true
        global _chol = chol
    catch e
        if isa(e, PosDefException)
            global _isPositiveDefinite = false
            rethrow(e)
        else
           rethrow(e)
        end
    end
    return
end

function _is_success()
    return _isOK
end

function _get_cond2()
    mat = _mat
    mat_norm = norm(mat)
    inv_norm = norm(inv(mat))
    return mat_norm * inv_norm
end
