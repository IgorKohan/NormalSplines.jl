module NormalSplines

#### Inteface deinition
export interpolate, evaluate
export RK_W1, RK_W2, RK_W3, RK_H1, RK_H2, RK_H3
export get_cond, get_epsilon
####

using LinearAlgebra

abstract type ReproducingKernel end
abstract type ReproducingKernel_1 <: ReproducingKernel   end
abstract type ReproducingKernel_2 <: ReproducingKernel_1 end
abstract type ReproducingKernel_3 <: ReproducingKernel_2 end

include("./ReproducingKernels.jl")
include("./GramMatrix.jl")
include("./Interpolate.jl")

"""
`interpolate(knots::Vector{T}, values::Vector{T}, kernel :: RK, factorize = true) where {T <: AbstractFloat, RK <: ReproducingKernel_1}`

Create an interpolating normal spline by values of function in knots.

# Arguments
- `knots::Vector{T}`: locations of given function values.
- `values::Vector{T}`: function values in `knots`.
- `kernel::RK`: reproducing kernel of Hilbert space the normal spline is constructed in.
                It must be a struct object of the following type:
                  `RK_W1` or `RK_H1` if the spline is constructing as a continuous function,
                  `RK_W2` or `RK_H2` if the spline is constructing as a differentiable function,
                  `RK_W3` or `RK_H3` if the spline is constructing as a twice differentiable function.
- `factorize::Bool`: defines if it is necessary to compute the Gram matrix factorization.
                  Must be set to `true` if the spline is constructing first time.
                  Can be set to `false` if the spline is constructing with the same knots
                  as it was done first time but with different function values in knots.

Returns : nothing.
"""
function interpolate(knots::Vector{T},
                     values::Vector{T},
                     kernel :: RK,
                     factorize = true) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
    if factorize
        _prepare_spline(kernel, knots)
    end
    _construct_spline(values)
    return
end

"""
`interpolate(knots::Vector{T}, values::Vector{T},
             d_knots::Vector{T}, d_values::Vector{T},
             kernel :: RK, derivative = 1, factorize = true) where {T <: AbstractFloat, RK <: ReproducingKernel_2}`

Create an interpolating normal spline by values of function
and it's first or second derivatives in knots.

# Arguments
- `knots::Vector{T}`: locations of given function values.
- `values::Vector{T}`: function values in `knots`.
- `d_knots::Vector{T}`: locations of given function first derivative values (if `derivative`=`1`)
                        or second derivative values (if `derivative`=`2`).
- `d_values::Vector{T}`: values of function first derivative (if derivative=`1`) or
                         values of function second derivative (if derivative=`2`) in `d_knots`.
- `kernel::RK`: reproducing kernel of Hilbert space the normal spline is constructed in.
                In a case when `derivative`=`1`, `kernel` must be a struct object of the following type:
                     `RK_W2` or `RK_H2` if the spline is constructing as a differentiable function,
                     `RK_W3` or `RK_H3` if the spline is constructing as a twice differentiable function.
                In a case when `derivative`=`2`, `kernel` must be a struct object of the following type:
                     `RK_W3` or `RK_H3` - the spline is constructing as a twice differentiable function,
- `derivative::Int`: `1` - `d_knots` and `d_values` provide data of function first derivative,
                     `2` - `d_knots` and `d_values` provide data of function second derivative.
- `factorize::Bool`: defines if it is necessary to compute the Gram matrix factorization.
                     Must be set to `true` if the spline is constructing first time.
                     Can be set to `false` if the spline is constructing with the same knots
                     as it was done first time but with different function values in knots.
Returns : nothing.
"""
function interpolate(knots::Vector{T},
                     values::Vector{T},
                     d_knots::Vector{T},
                     d_values::Vector{T},
                     kernel :: RK,
                     derivative = 1,
                     factorize = true) where {T <: AbstractFloat, RK <: ReproducingKernel_2}
    if derivative == 1
        if factorize
            _prepare_spline_d1(kernel, knots, d_knots)
        end
        _construct_spline_d1(values, d_values)
    elseif derivative == 2
        if isa(kernel, RK_H2) || isa(kernel, RK_W2)
           Throw(ArgumentError("kernel: Incorrect parameter type (must be RK_H3 or RK_W3)."))
        end
        if factorize
            _prepare_spline_d2(kernel, knots, d_knots)
        end
        _construct_spline_d2(values, d_values)
    else
        Throw(DomainError(derivative, "'derivative' parameter value must be 1 or 2."))
    end

    return
end

"""
`interpolate(knots::Vector{T}, values::Vector{T},
             d1_knots::Vector{T}, d1_values::Vector{T},
             d2_knots::Vector{T}, d2_values::Vector{T},
             kernel :: RK, factorize = true) where {T <: AbstractFloat, RK <: ReproducingKernel_3}`

Create an interpolating normal spline by values of function
and it's first and second derivatives in knots.

# Arguments
- `knots::Vector{T}`: locations of given function values.
- `values::Vector{T}`: function values in `knots`.
- `d1_knots::Vector{T}`: locations of given function derivative values.
- `d1_values::Vector{T}`: values of function derivative in `d1_knots`.
- `d2_knots::Vector{T}`: locations of given function second derivative values.
- `d2_values::Vector{T}`: values of function second derivative in `d2_knots`.
- `kernel::RK`: reproducing kernel of Hilbert space the normal spline is constructed in.
                It must be a struct object of the following type:
                  `RK_W3` or `RK_H3` - the spline is constructing as a twice differentiable function.
- `factorize::Bool`: defines if it is necessary to compute the Gram matrix factorization.
                  Must be set to `true` if the spline is constructing first time.
                  Can be set to `false` if the spline is constructing with the same knots
                  as it was done first time but with different function values in knots.

Returns : nothing.
"""
function interpolate(knots::Vector{T},
                     values::Vector{T},
                     d1_knots::Vector{T},
                     d1_values::Vector{T},
                     d2_knots::Vector{T},
                     d2_values::Vector{T},
                     kernel :: RK,
                     factorize = true) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
    if factorize
        _prepare_spline_d12(kernel, knots, d1_knots, d2_knots)
    end
    _construct_spline_d12(values, d1_values, d2_values)
    return
end

"""
`evaluate(points::Vector{T}, derivative::Int = 0) where T <: AbstractFloat`

Evaluate normal spline and it's derivatives.

# Arguments
- `points::Vector{T}`: locations at which spline or its derivative are evaluating.
- `derivative::Int`: `0` - spline evaluation.
                     `1` - first derivative of spline evaluation.
                     `2` - second derivative of spline evaluation.

Returns: `Vector{T}` of spline values or its derivatives at the locations
          defined in `points`.
"""
function evaluate(points::Vector{T}, derivative::Int = 0) where T <: AbstractFloat
    if derivative == 0
        return _evaluate(points)
    elseif derivative == 1
        return _evaluate_d1(points)
    elseif derivative == 2
        return _evaluate_d2(points)
    else
        Throw(DomainError(derivative, "'derivative' parameter value must be 0, 1 or 2."))
    end
end

"""
`evaluate(point::T, derivative::Int = 0) where T <: AbstractFloat`

Evaluate normal spline and it's derivatives.

# Arguments
- `point::T`: location at which spline or its derivative are evaluating.
- `derivative::Int`: `0` - spline evaluation.
                     `1` - first derivative of spline evaluation.
                     `2` - second derivative of spline evaluation.

Returns: spline value or its derivative at the `point` location.
"""
function evaluate(point::T, derivative::Int = 0) where T <: AbstractFloat
    v_points = Vector{T}(undef, 1)
    v_points[1] = point
    return evaluate(v_points, derivative)[1]
end

"""
`get_cond()`

Estimate condition number of Gram matrix.
"""
function get_cond()
    return _get_cond()
end

"""
`get_epsilon()`

Get value of parameter `ε` of struct object `RK_H1`,
`RK_H2` or `RK_H3`.
"""
function get_epsilon()
    return _get_epsilon()
end

end # module
