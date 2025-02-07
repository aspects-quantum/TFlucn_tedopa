

 using DrWatson, Plots, LaTeXStrings

gr()  # Use GR backend, but you can switch to PyPlot if needed
function ticks_length!(;tl=0.02)
    p = Plots.current()
    xticks, yticks = Plots.xticks(p)[1][1], Plots.yticks(p)[1][1]
    xl, yl = Plots.xlims(p), Plots.ylims(p)
    x1, y1 = zero(yticks) .+ xl[1], zero(xticks) .+ yl[1]
    sz = p.attr[:size]
    r = sz[1]/sz[2]
    dx, dy = tl*(xl[2] - xl[1]), tl*r*(yl[2] - yl[1])
    plot!([xticks xticks]', [y1 y1 .+ dy]', c=:black, labels=false)
    plot!([x1 x1 .+ dx]', [yticks yticks]', c=:black, labels=false, xlims=xl, ylims=yl)
    return Plots.current()
end
# Define data
a_plus = [0.09495992168234713, 0.33702123217274094, 0.6372444896736239, 0.9231367705835735,
    1.1615464001371252, 1.347733484238672, 1.4890100635883856, 1.595849570599494,
    1.677332024145704, 1.7402102037323495, 1.790046901376592, 1.8304485200765086,
    1.8634051422401667, 1.890939939037258, 1.9149037738675725, 1.93495135016066,
    1.9527759838347123, 1.9679265433206998, 1.9808005761639165, 1.9926380385688867,
    2.002799904210332, 2.012851674669058, 2.0217916794514443, 2.030263953987257,
    2.038063567343079, 2.044545032604059, 2.0505445534867364, 2.0575423103798114,
    2.0640886354658194, 2.0696213998269517, 2.0751770882007134]

a_up = [0.09552011750390858, 0.34278645428748716, 0.6586830462975073, 0.9725049773940453,
    1.2483702020746037, 1.477447258343144, 1.6632476355539043, 1.813411162889399,
    1.9355827522857658, 2.035944679745456, 2.1189882271626956, 2.188229279321799,
    2.2468506799215917, 2.2974034397288943, 2.341680621541797, 2.3795126059699805,
    2.4135375490451954, 2.4439780949331515, 2.4708301922974116, 2.4948263635099726,
    2.5164834344953952, 2.5357876912511426, 2.553564175353163, 2.5699792516313194,
    2.5851434156168316, 2.5997227095899067, 2.613347278374249, 2.6258577431914243,
    2.637438606399008, 2.648724372118132, 2.658218253384091]

a_in = [0.024466749688153988, 0.20500971482632774, 0.49715656024986365, 0.8176927994838703,
    1.1133439355639254, 1.3636660084147194, 1.5663140181582313, 1.7272349629271395,
    1.8547726635280228, 1.9561300224757139, 2.0373350229850016, 2.1032705976393467,
    2.1571365328331393, 2.2015699058764757, 2.23870318478782, 2.269818758464505,
    2.2960795674929733, 2.3186733585261092, 2.3379718925965296, 2.354656820179843,
    2.3692663035069925, 2.3822852229335547, 2.3937976873363755, 2.404002163018495,
    2.4131111441075084, 2.421474585636393, 2.429016870453073, 2.435708932073047,
    2.441760765034304, 2.4472897125763287, 2.4522404292131528]

# Time step duration
tau = 0.04
time_steps = collect(0:length(a_plus)-1) * tau

# Plot with adjusted settings
p1 = plot()

plot!(time_steps, real(a_plus), label = L"\ \ 0.5   \quad \quad  1   \quad \quad   0   \quad \quad   |+\rangle\langle +|",
    seriesalpha = 0.5, linewidth = 8, legendfontsize = 12, size = (900, 700))

plot!(time_steps, real(a_in), label = L"\ \  0.5  \quad \quad   0   \quad \quad   1  \quad \quad  |\uparrow\rangle\langle \uparrow |",
    seriesalpha = 0.5, linewidth = 8)

plot!(time_steps, real(a_up), label = L"\ \ 0.5   \quad \quad  1   \quad \quad   0   \quad \quad  |\uparrow\rangle\langle \uparrow |",
    seriesalpha = 0.5, linewidth = 8)

# Customize axes labels
xlabel!(L"t", fontsize = 20)
ylabel!(L"\langle Q \rangle", fontsize = 20)

# Set tick positions and convert tick labels to LaTeX strings automatically
xticks = range(0, stop = maximum(time_steps), length = 5)
yticks = range(0, stop = maximum(a_plus) + 0.2, length = 5)

xtick_labels = [L"\$" * string(round(x, digits = 2)) * "\$" for x in xticks]
ytick_labels = [L"\$" * string(round(y, digits = 2)) * "\$" for y in yticks]

# Set tick and label font sizes, and apply LaTeX tick labels
plot!(xticks = (xticks, xtick_labels), yticks = (yticks, ytick_labels),
    xtickfont = font(30), ytickfont = font(30), guidefont = font(30))

ticks_length!(tl=.04)

# Customize plot appearance, place legend inside
plot!(legend = :right, grid = false, framestyle = :box, legendposition = (.4, 0.5))

# Add a legend title with parameter names, acting as column headers
legend_title = L"\quad \alpha \quad \quad \Omega \quad  \quad   \omega_0   \quad \quad   \rho_{S}(0)"
plot!(legendtitle = legend_title, gridlinewidth = 2,  legendfontsize = 25, legendtitlefontsize = 25)


# Display and save the plot
display(p1)
file_name_png = "improved_plot_legend.png"
safesave(file_name_png, p1)
