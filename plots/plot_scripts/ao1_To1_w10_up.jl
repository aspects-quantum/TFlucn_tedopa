

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


#T = 0.1, alpha = 0.1, N_chain = 140, maxdim = 80, cutoff = -12, tau = 0.003, jump = 10, boson_dim = 12, omega = 10;
mQ1 = [-4.322974185804232e-10;
0.08259796187253249;
0.26799643188509636;
0.4610454825483667;
0.576925797500231;
0.6749189431007785;
0.7962810236997417;
0.8846211159587202;
0.9415341057886427;
0.9855039715505398;
1.0108764601683118;
1.0242244663611586;
1.0484870829892419;
1.0673332577582153;
1.0724122084383296;
1.085596845043796;
1.103953926312762;
1.1106237702551813;
1.1200633733181404;
1.1356791812989122;
1.146429523794396;
1.1346463718706181;
1.171108797984156;
1.1983172295769715;
1.194336980520428;
1.166796268005421;
1.174559781772046;
1.200943683395714;
1.2079340687446296;
1.2098608979149998;
1.2477195821541123;
1.2622698549314157;

]
vQ1 = [2.031575586072204e-19;
2.430377527300094;
7.498635533842515;
11.557487014279792;
11.638818170016382;
11.767166971926478;
10.776695067290532;
12.741918695930112;
13.697074019770278;
13.449796682417048;
13.412077323303393;
12.762480496924892;
12.951944463110447;
12.582703495576608;
12.274228067313501;
12.057911462254369;
11.975117572174048;
11.157368355852363;
10.75909089532909;
10.735025335905343;
10.985281524464485;
10.218090209924775;
9.757081673375732;
9.704889107827098;
7.686269202010412;
6.217786526037488;
8.946894964274437;
10.510296130878343;
9.683775408865996;
9.88678339549587;
13.968126439512266;
16.829141613121003;
]
#

#= 

mQ1 = mQ1[1:30]
vQ1 = vQ1[1:30]
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
