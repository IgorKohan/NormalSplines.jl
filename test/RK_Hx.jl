x = [0.0, 1.0, 2.0]                        # function knots
u = x.^2                                   # function values in knots
s = [2.0]                                  # function first derivative knots
v = [4.0]                                  # function first derivative values
t = [0.0, 1.0]                             # function second derivative knots
w = [2.0 ,2.0]                             # function second derivative values
d1_u = v                                   # exact function first derivative values in knots
d2_u = w                                   # exact function second derivative values in knots

@testset "Bessel potentials space reproducing kernels" begin
    @testset "RK_H1 kernel" begin
        interpolate(x, u, RK_H1(0.1))  # create spline
        cond = get_cond()              # get estimation of the Gram matrix condition number
        @test cond ≈ 100.0

        σ = evaluate(x)                # evaluate spline in knots
        @test σ ≈ u                    # compare with exact function values in knots

        # Check that we get close when evaluating near the knots
        p = x .+ 1e-3*randn(size(x))   # evaluation points near the knots
        f = p.^2                       # exact function values in evaluation points
        σ = evaluate(p)                # evaluate spline in evaluation points
        # compare spline values with exact function values in evaluation point
        @test all(isapprox.(σ, f, atol = 0.05))
    end

    @testset "RK_H2 kernel" begin
        interpolate(x, u, RK_H2(0.1))  # create spline
        cond = get_cond()              # get estimation of the gram matrix condition number
        @test cond ≈ 1.0e5

        σ = evaluate(x)                # evaluate spline in knots
        @test σ ≈ u                    # compare with exact function values in knots

        # Check that we get close when evaluating near the knots
        p = x .+ 1e-3*randn(size(x))   # evaluation points near the knots
        f = p.^2                       # exact function values in evaluation points
        σ = evaluate(p)                # evaluate spline in evaluation points
        # compare spline values with exact function values in evaluation point
        @test all(isapprox.(σ, f, atol = 1e-2))

        ###
        interpolate(x, u, s, v, RK_H2(0.1)) # create spline by function and
                                            # first derivative values in knots
        cond = get_cond()              # get estimation of the gram matrix condition number
        @test cond ≈ 1.0e5

        σ = evaluate(x)                # evaluate spline in knots
        @test σ ≈ u                    # compare with exact function values in knots

        # Check that we get close when evaluating near the knots
        p = x .+ 1e-3*randn(size(x))   # evaluation points near the knots
        f = p.^2                       # exact function values in evaluation points
        σ = evaluate(p)                # evaluate spline in evaluation points
        # compare spline values with exact function values in evaluation point
        @test all(isapprox.(σ, f, atol = 1e-2))

        d1_σ = evaluate(s, 1)             # evaluate spline first derivative in knots
        @test d1_σ ≈ d1_u                 # compare with exact function first derivative values in knots

        # Check that we get close when evaluating near the knots
        ps = s .+ 1e-3*randn(size(s))   # evaluation points near the knots
        d1_f = ps.*2                    # exact function first derivative values in evaluation points
        d1_σ = evaluate(ps, 1)             # evaluate spline first derivative in evaluation points
        # compare spline first derivative values with exact function first derivative values in evaluation point
        @test all(isapprox.(d1_σ, d1_f, atol = 1e-2))
    end

    @testset "RK_H3 kernel" begin
        interpolate(x, u, RK_H3(0.1))  # create spline
        cond = get_cond()              # get estimation of the gram matrix condition number
        @test cond ≈ 1.0e7

        σ = evaluate(x)                # evaluate spline in knots
        @test σ ≈ u                    # compare with exact function values in knots

        # Check that we get close when evaluating near the knots
        p = x .+ 1e-3*randn(size(x))   # evaluation points near the knots
        f = p.^2                       # exact function values in evaluation points
        σ = evaluate(p)                # evaluate spline in evaluation points
        # compare spline values with exact function values in evaluation point
        @test all(isapprox.(σ, f, atol = 1e-2))

        ###
        interpolate(x, u, s, v, RK_H3(0.1))     # create spline
        cond = get_cond()              # get estimation of the gram matrix condition number
        @test cond ≈ 1.0e8

        σ = evaluate(x)                # evaluate spline in knots
        @test σ ≈ u                    # compare with exact function values in knots

        # Check that we get close when evaluating near the knots
        p = x .+ 1e-3*randn(size(x))   # evaluation points near the knots
        f = p.^2                       # exact function values in evaluation points
        σ = evaluate(p)                # evaluate spline in evaluation points
        # compare spline values with exact function values in evaluation point
        @test all(isapprox.(σ, f, atol = 1e-2))

        d1_σ = evaluate(s, 1)             # evaluate spline first derivative in knots
        @test d1_σ ≈ d1_u                 # compare with exact function first derivative values in knots

        # Check that we get close when evaluating near the knots
        ps = s .+ 1e-3*randn(size(s))   # evaluation points near the knots
        d1_f = ps.*2                    # exact function first derivative values in evaluation points
        d1_σ = evaluate(ps, 1)             # evaluate spline first derivative in evaluation points
        # compare spline first derivative values with exact function first derivative values in evaluation point
        @test all(isapprox.(d1_σ, d1_f, atol = 1e-2))

        ###
        interpolate(x, u, t, w, RK_H3(0.1), 2)     # create spline
        cond = get_cond()              # get estimation of the gram matrix condition number
        @test cond ≈ 1.0e8

        σ = evaluate(x)                # evaluate spline in knots
        @test σ ≈ u                    # compare with exact function values in knots

        # Check that we get close when evaluating near the knots
        p = x .+ 1e-3*randn(size(x))   # evaluation points near the knots
        f = p.^2                       # exact function values in evaluation points
        σ = evaluate(p)                # evaluate spline in evaluation points
        # compare spline values with exact function values in evaluation point
        @test all(isapprox.(σ, f, atol = 1e-2))

        d1_σ = evaluate(s, 1)          # evaluate spline first derivative in knots
        @test d1_σ ≈ d1_u atol = 0.05  # compare with exact function first derivative values in knots

        # Check that we get close when evaluating near the knots
        ps = s .+ 1e-3*randn(size(s))   # evaluation points near the knots
        d1_f = ps.*2                    # exact function first derivative values in evaluation points
        d1_σ = evaluate(ps, 1)          # evaluate spline first derivative in evaluation points
        # compare spline first derivative values with exact function first derivative values in evaluation point
        @test all(isapprox.(d1_σ, d1_f, atol = 0.05))

        d2_σ = evaluate(t, 2)           # evaluate spline second derivative in knots
        @test d2_σ ≈ d2_u               # compare with exact function second derivative values in knots

        # Check that we get close when evaluating near the knots
        pt = t .+ 1e-3*randn(size(t))   # evaluation points near the knots
        d2_f = 2.0                      # exact function second derivative values in evaluation points
        d2_σ = evaluate(pt, 2)          # evaluate spline second derivative in evaluation points
        # compare spline second derivative values with exact function second derivative values in evaluation point
        @test all(isapprox.(d2_σ, d2_f, atol = 1e-2))

        ###
        interpolate(x, u, s, v, t, w, RK_H3(0.1)) # create spline
        cond = get_cond()              # get estimation of the gram matrix condition number
        @test cond ≈ 1.0e9

        σ = evaluate(x)                # evaluate spline in knots
        @test σ ≈ u                    # compare with exact function values in knots

        # Check that we get close when evaluating near the knots
        p = x .+ 1e-3*randn(size(x))   # evaluation points near the knots
        f = p.^2                       # exact function values in evaluation points
        σ = evaluate(p)                # evaluate spline in evaluation points
        # compare spline values with exact function values in evaluation point
        @test all(isapprox.(σ, f, atol = 1e-2))

        d1_σ = evaluate(s, 1)          # evaluate spline first derivative in knots
        @test d1_σ ≈ d1_u              # compare with exact function first derivative values in knots

        # Check that we get close when evaluating near the knots
        ps = s .+ 1e-3*randn(size(s))  # evaluation points near the knots
        d1_f = ps.*2                   # exact function first derivative values in evaluation points
        d1_σ = evaluate(ps, 1)         # evaluate spline first derivative in evaluation points
        # compare spline first derivative values with exact function first derivative values in evaluation point
        @test all(isapprox.(d1_σ, d1_f, atol = 1e-2))

        d2_σ = evaluate(t, 2)          # evaluate spline second derivative in knots
        @test d2_σ ≈ d2_u              # compare with exact function second derivative values in knots

        # Check that we get close when evaluating near the knots
        pt = t .+ 1e-3*randn(size(t))   # evaluation points near the knots
        d2_f = 2.0                      # exact function second derivative values in evaluation points
        d2_σ = evaluate(pt, 2)          # evaluate spline second derivative in evaluation points
        # compare spline second derivative values with exact function second derivative values in evaluation point
        @test all(isapprox.(d2_σ, d2_f, atol = 1e-2))
    end

end
