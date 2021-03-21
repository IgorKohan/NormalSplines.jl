export demo

using Gadfly

function demo()
# Setup
    x = collect(1.0:1.0:20)
    u = x.*0.0
    s = x
    v = x.*0.0
    t = x
    w = x.*0.0

    for i in 6:10
        u[i] = 1.0
    end
    for i in 11:14
        u[i] = -0.2 * i + 3.0
    end
    for i in 11:14
        v[i] = -0.2
    end

    p = collect(1.0:0.2:20)
    r = p.*0.0
    tol = 1e-7
    i = 0
    for pi in p
        i += 1
        if pi >= (6.0 - tol) && pi <= (10.0 + tol)
            r[i] = 1.0
        end
        if pi > (10.0 + tol) && pi < (15.0 - tol)
            r[i] = -0.2 * p[i] + 3.0
        end
    end
####
    interpolate(x, u, RK_W3())

    cond = get_cond()
    ε = get_epsilon()
    println("cond = $cond, ε = $ε")

    σ = evaluate(p)

    δ = σ - r
    rms = round(norm(δ)/sqrt(length(δ)); digits=2)
    println("rms = $rms")

    set_default_plot_size(13cm, 13cm)
    plt = plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
          layer(x = p, y = r, Geom.line, Theme(default_color=colorant"red")),
          layer(x = p, y = σ, Geom.line, Theme(default_color=colorant"blue")),
          Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
          Guide.manual_color_key("Legend", ["Points", "True", "Spline"], ["orange", "red", "blue"]),
          Guide.xlabel("knots, points"), Guide.ylabel("f, σ"),
          Guide.title("Fig.1a"))
    draw(SVG("plot_W3.svg", 13cm, 13cm), plt)

    σ = evaluate(p, 1)
    plt = plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
          layer(x = p, y = σ, Geom.line, Theme(default_color=colorant"blue")),
          Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
          Guide.manual_color_key("Legend", ["Points", "Spline 1st derivative"], ["orange", "blue"]),
          Guide.xlabel("knots, points"), Guide.ylabel(raw"σ'"),
          Guide.title("Fig.2a"))
    draw(SVG("plot_W3_d1.svg", 13cm, 13cm), plt)

    σ = evaluate(p, 2)
    plt = plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
          layer(x = p, y = σ, Geom.line, Theme(default_color=colorant"blue")),
          Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
          Guide.manual_color_key("Legend", ["Points", "Spline 2nd derivative"], ["orange", "blue"]),
          Guide.xlabel("knots, points"), Guide.ylabel(raw"σ''"),
          Guide.title("Fig.3a"))
    draw(SVG("plot_W3_d2.svg", 13cm, 13cm), plt)

###
      interpolate(x, u, s, v, t, w, RK_W3())

      cond = get_cond()
      ε = get_epsilon()
      println("cond = $cond, ε = $ε")

      σ = evaluate(p)

      δ = σ - r
      rms = round(norm(δ)/sqrt(length(δ)); digits=2)
      println("rms = $rms")

      set_default_plot_size(13cm, 13cm)
      plt = plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
            layer(x = p, y = r, Geom.line, Theme(default_color=colorant"red")),
            layer(x = p, y = σ, Geom.line, Theme(default_color=colorant"blue")),
            Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
            Guide.manual_color_key("Legend", ["Points", "True", "Spline"], ["orange", "red", "blue"]),
            Guide.xlabel("knots, points"), Guide.ylabel("f, σ"),
            Guide.title("Fig.1b"))
      draw(SVG("plot_d12_W3.svg", 13cm, 13cm), plt)

      σ = evaluate(p, 1)
      plt = plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
            layer(x = p, y = σ, Geom.line, Theme(default_color=colorant"blue")),
            Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
            Guide.manual_color_key("Legend", ["Points", "Spline 1st derivative"], ["orange", "blue"]),
            Guide.xlabel("knots, points"), Guide.ylabel(raw"σ'"),
            Guide.title("Fig.2b"))
      draw(SVG("plot_d12_W3_d1.svg", 13cm, 13cm), plt)

      σ = evaluate(p, 2)
      plt = plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
            layer(x = p, y = σ, Geom.line, Theme(default_color=colorant"blue")),
            Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
            Guide.manual_color_key("Legend", ["Points", "Spline 2nd derivative"], ["orange", "blue"]),
            Guide.xlabel("knots, points"), Guide.ylabel(raw"σ''"),
            Guide.title("Fig.3b"))
      draw(SVG("plot_d12_W3_d2.svg", 13cm, 13cm), plt)

end
