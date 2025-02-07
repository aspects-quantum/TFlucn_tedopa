

using DrWatson, Plots, LaTeXStrings

gr()  # Use GR backend, but you can switch to PyPlot if needed
function ticks_length!(; tl = 0.02)
	p = Plots.current()
	xticks, yticks = Plots.xticks(p)[1][1], Plots.yticks(p)[1][1]
	xl, yl = Plots.xlims(p), Plots.ylims(p)
	x1, y1 = zero(yticks) .+ xl[1], zero(xticks) .+ yl[1]
	sz = p.attr[:size]
	r = sz[1] / sz[2]
	dx, dy = tl * (xl[2] - xl[1]), tl * r * (yl[2] - yl[1])
	plot!([xticks xticks]', [y1 y1 .+ dy]', c = :black, labels = false)
	plot!([x1 x1 .+ dx]', [yticks yticks]', c = :black, labels = false, xlims = xl, ylims = yl)
	return Plots.current()
end
# Define data


#T = 0.01, alpha = 0.5, N_chain = 160, maxdim = 80, cutoff = -12, tau = 0.005, boson_dim = 12;

mQ1 = [0; 0.03841308361292526;
0.146843005811204;
0.30776743850586313;
0.4992163119644793;
0.7009292450111113;
0.8983543679879833;
1.0821796910641488;
1.2478154379677324;
1.3943425265324039;
1.5217714010211976;
1.6320243166349606;
1.7279303345868304;
1.811289360590564;
1.8834458527589595;
1.9460981013886252;
2.0014873602357035;
2.049342302311088;
2.090475510788706;
2.125903565949205;
2.1595039412065336;
2.187687782315866;
2.210128979172345;
2.230231869697804;
2.24772634597866;
]
vQ1 = [0; 0.5683140067308483;
2.1216902624195084;
4.319419044185123;
6.730185043111731;
9.039532476568931;
11.056853268355265;
12.719427722282338;
14.062821551866461;
15.10264829740894;
15.916635760770642;
16.542006659667152;
17.04023777479824;
17.3652386034426;
17.637239333428184;
17.807929429178678;
17.9954303209052;
18.14509269370415;
18.241550225795635;
18.32156667682286;
18.33818216308376;
18.3906199246966;
18.363109587626884;
18.234909781506353;
18.314401833791706;
]
#

#= 

mQ1 = mQ1[1:43]
vQ1 = vQ1[1:43]
 =#
# Time step duration
tau = 0.005
nt = 250
ttotal = nt * tau  # Total time evolution
#time_steps = collect(0:10*tau:ttotal)

time_steps = collect(0:lastindex(mQ1)-1) * 5 * tau

# Set tick positions and convert tick labels to LaTeX strings automatically
xticks = range(0, stop = maximum(time_steps), length = 5)
yticks1 = range(0, stop = 6, length = 5)
yticks2 = range(0, stop = 60, length = 5)
#= 
xtick_labels = [L"\$" * string(round(x, digits = 2)) * "\$" for x in xticks]
ytick_labels1 = [L"\$" * string(round(y, digits = 1)) * "\$" for y in yticks1]
ytick_labels2 = [L"\$" * string(round(y, digits = 1)) * "\$" for y in yticks2] =#

xtick_labels = [string(round(x, digits = 2)) for x in xticks]
ytick_labels1 = [string(round(y, digits = 1)) for y in yticks1]
ytick_labels2 = [string(round(y, digits = 1)) for y in yticks2]


p1 = plot(time_steps, mQ1, xaxis = L"t", yaxis = L"\langle Q \rangle", 
xtickfont = font(12), ytickfont = font(12), xticks = (xticks, xtick_labels),# yticks = (yticks1, ytick_labels1), 
xguidefontsize=18,yguidefontsize=20, color = :orange,
seriesalpha = 0.9, linewidth = 6, label = L"\ \alpha = 0.5")


plot!(twinx(),time_steps, vQ1, line = :dot, color = :orange, ytickfont = font(12),
seriesalpha = 0.9, linewidth = 6, label = false, yguidefontsize=20, yaxis = L"\langle\langle Q^2 \rangle\rangle")

plot!(widen = false)
vline!([xlims(p1)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p1)[2]], lc = :black, lw = 2, label = false)
#plot!(legend = :outertopright)

title!(L"\epsilon = 1, \Delta = 0, T = 0.01")

#ticks_length!(tl=.04)

# Customize plot appearance, place legend inside
plot!(legend = :right, grid = false, legendfontsize = 16)

display("image/png", p1)
