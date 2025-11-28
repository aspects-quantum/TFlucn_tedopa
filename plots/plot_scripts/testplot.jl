using DrWatson          # fine but not needed here
using Plots
using LaTeXStrings
using Plots.PlotMeasures
gr()                    # or your preferred backend

function broken_y_plot(x, y_low, y_high;
                       ylim_low  = (0, 2),
                       ylim_high = (38, 42),
                       xlabel    = L"x",
                       ylabel    = L"y",
                       title     = "")

    # --- upper segment (large values) ---
    p_top = scatter(
        x, y_high;
        ylim = ylim_high,
        legend = false,
        xaxis = false,      # hide x-axis and ticks
        grid = false,
        framestyle = :box,
        ylabel     = ylabel,
        title      = title,
        bottom_margin = -5mm,    # pull panels closer
        markerstrokewidth = 0,
    )

    # --- lower segment (small values) ---
    p_bot = scatter(
        x, y_low;
        ylim       = ylim_low,
        legend     = false,
        #framestyle = :box,
        grid = false,
        xlabel     = xlabel,
        ylabel     = ylabel,
        top_margin = -5mm,
        markerstrokewidth = 0,
    )

    # combine so it *looks* like a single plot with a y-break
    plot(p_top, p_bot;
         layout = @layout([a; b]),
         link   = :x,
         size   = (500, 500),
         margin = 3mm)
end

x      = 1:6
y_low  = [0.3, 0.8, 1.2, 0.6, 1.5, 0.9]          # near 0–2
y_high = [40.0, 39.5, 41.0, 40.5, 39.8, 40.2]    # way up at ~40


plt = broken_y_plot(
    x, y_low, y_high;
    ylim_low  = (0, 2),
    ylim_high = (38, 42),
    xlabel    = L"t",
    ylabel    = L"Q",
    title     = "Broken y-axis demo"
)

display(plt)

