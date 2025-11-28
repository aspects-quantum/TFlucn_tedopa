using Plots
using LaTeXStrings
using Plots.PlotMeasures
gr()

"""
    broken_y_plot(x, y_low, y_high;
                  ylim_low=(0,2), ylim_high=(38,42),
                  xlabel=L"x", ylabel=L"y", title="")

Two vertical segments with a broken y-axis:
- lower panel shows `ylim_low`
- upper panel shows `ylim_high`
- only bottom x-axis is visible
- a small 'wiggle' is drawn at the break
"""
function broken_y_plot(x, y_low, y_high;
                       ylim_low  = (0, 2),
                       ylim_high = (38, 42),
                       xlabel    = L"x",
                       ylabel    = L"y",
                       title     = "")

    # --- upper segment (large values) ---
    p_top = scatter(
        x, y_high;
        ylim       = ylim_high,
        legend     = false,
        framestyle = :left,
        xaxis      = false,              # hide x axis completely
        xticks     = false,
        xguide     = "",
        #xforeground_color = :white,
        ylabel     = ylabel,
        title      = title,
        bottom_margin = -5mm,
        markerstrokewidth = 0,
    )

    # --- lower segment (small values) ---
    p_bot = scatter(
        x, y_low;
        ylim       = ylim_low,
        legend     = false,
        #framestyle = :top,
        xlabel     = xlabel,
        ylabel     = ylabel,
        top_margin = -5mm,
        markerstrokewidth = 0,
    )

    # ===== add the little wiggly break marks =====
    # we'll draw a small 'X' / zig-zag on the left side of each panel
    # position in x (slightly to the left of the first data point)
    x0 = minimum(x) - 0.1
    # small vertical sizes, relative to each panel's y-range
    δy_top = 0.02 * (ylim_high[2] - ylim_high[1])
    δy_bot = 0.02 * (ylim_low[2]  - ylim_low[1])

    # top panel: wiggle at bottom of its range (break to lower panel)
    y_top = ylim_high[1]
    plot!(p_top,
          [x0-0.03, x0+0.03], [y_top+δy_top, y_top-δy_top],
          lw=2, color=:black)
    plot!(p_top,
          [x0-0.03, x0+0.03], [y_top-δy_top, y_top+δy_top],
          lw=2, color=:black)

    # bottom panel: wiggle at top of its range (break to upper panel)
    y_bot = ylim_low[2]
    plot!(p_bot,
          [x0-0.03, x0+0.03], [y_bot+δy_bot, y_bot-δy_bot],
          lw=2, color=:black)
    plot!(p_bot,
          [x0-0.03, x0+0.03], [y_bot-δy_bot, y_bot+δy_bot],
          lw=2, color=:black)

    # combine panels
    plot(p_top, p_bot;
         layout = @layout([a; b]),
         link   = :x,
         size   = (500, 500),
         margin = 3mm)
end
x      = 1:6
y_low  = [0.3, 0.8, 1.2, 0.6, 1.5, 0.9]
y_high = [40.0, 39.5, 41.0, 40.5, 39.8, 40.2]

plt = broken_y_plot(
    x, y_low, y_high;
    ylim_low  = (0, 2),
    ylim_high = (38, 42),
    xlabel    = L"t",
    ylabel    = L"Q",
    title     = "Broken y-axis demo"
)
display(plt)
