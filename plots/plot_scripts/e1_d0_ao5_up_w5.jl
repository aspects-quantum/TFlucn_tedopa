

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

mQ1 = [0;
0.024780033331352848;
0.09653322605835667;
0.2081729246322872;
0.348488747324072;
0.5086833768285567;
0.6787604909099137;
0.8506584149375362;
1.0188813075094045;
1.1369961854096375;
1.2462493541743545;
1.343859944979723;
1.553258106553309;
1.7221886754519047;
1.8473222160597618;
1.9391445713186708;
2.0182279447785683;
2.093656989447123;
2.1679788528461295;
2.238272194150845;
2.320336883706607;
2.3700892902292727;
2.4044355839462597;
2.4387310930187343;
2.47530498220248;
2.5155905270363377;
2.561217530484358;
2.6019421881804923;
2.6192324487340364;
2.63327239441546;
2.646473010640129;
2.673072601233132;
2.7264501554548204;
2.7391148691876377;
2.748329220893723;
2.768161847778117;
2.7910701468094694;
2.801806016672613;
2.8205566960097816;

]
vQ1 = [0; 0.0006144964691742478;
0.0093254683177666;
0.043367242380106946;
5.161318066388364;
7.558611650274598;
9.787142890201375;
12.242982890452245;
13.413651968211898;
13.82263151150789;
13.180697965874518;
12.293173962405039;
15.58454122918483;
17.182459530481655;
18.197413051534223;
18.688977561197017;
18.708569865864874;
20.492379918254375;
20.71501604334855;
21.231176830972732;
21.431193702655563;
21.66024371508444;
21.93344060814882;
21.867724720623173;
21.41833156326507;
21.591410490123586;
22.199403355648947;
22.271477652026668;
21.747080995409526;
21.5235109627423;
21.21493327487693;
21.312512317838745;
22.075636245598762;
21.579447524409318;
21.014150832133925;
21.000567090441187;
21.095759502166278;
20.731071221555506;
20.622390957834593;
]
#



mQ1 = mQ1[1:39]
vQ1 = vQ1[1:39]


# Time step duration
tau = 0.002
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
seriesalpha = 0.9, linewidth = 6, label = L"\ \alpha = 0.5")


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
