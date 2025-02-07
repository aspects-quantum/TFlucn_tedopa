

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

mQ1 = [0; 0.11538137804578356;
0.44314407625446;
0.9348445828382738;
1.5287164820795502;
2.1653194874678094;
2.8006578157958546;
3.399849686288985;
3.9476874006968243;
4.43749703014376;
4.870338111923426;
5.246980090006293;
5.574183571294299;
5.850399231838691;
6.09493526000804;
6.269522152547097;
6.477660862970822;
6.627020597894135;
6.763840229472955;
6.870561606022656;
6.967775178926816;
7.053379485427154;
7.133471266665098;
7.185625180403613;
7.237964339885265;
7.279845107378538;
7.315794209646597;
7.3671283931431395;
7.405684097026345;
7.432727008987003;
]
vQ1 = [0; 1.7224293467574319;
6.783728965328398;
14.359022227258023;
22.577679413702825;
31.45832715832337;
40.054412349338584;
47.86643526717027;
54.372801677540345;
60.08870183908704;
65.00848215610566;
68.8936150271005;
72.52001013307034;
74.96632939579217;
77.34490324713961;
77.66068943500107;
81.36787568849843;
82.57836492833215;
84.59430910029371;
85.81938870759552;
86.86631289948006;
87.70916014639275;
88.48021793863751;
88.99337965550815;
89.50783234090098;
89.81805691083899;
90.11230937041267;
91.16838090972186;
91.90474624399715;
92.19156838707022;

]
#



mQ1 = mQ1[1:30]
vQ1 = vQ1[1:30]


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
xguidefontsize=18,yguidefontsize=20, color = :brown,
seriesalpha = 0.9, linewidth = 6, label = L"\ \alpha = 1.5")


plot!(twinx(),time_steps, vQ1, line = :dot, color = :brown, ytickfont = font(12),
seriesalpha = 0.9, linewidth = 6, label = false, yguidefontsize=20, yaxis = L"\langle\langle Q^2 \rangle\rangle")

plot!(widen = false)
vline!([xlims(p1)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p1)[2]], lc = :black, lw = 2, label = false)
#plot!(legend = :outertopright)

title!(L"\epsilon = 1, \Delta = 0, T = 0.01, \omega_C = 5, |1\rangle")

#ticks_length!(tl=.04)

# Customize plot appearance, place legend inside
plot!(legend = :right, grid = false, legendfontsize = 16)

display("image/png", p1)
