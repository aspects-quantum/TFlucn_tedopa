

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


#T = 0.1, alpha = 1.0, N_chain = 160, maxdim = 80, cutoff = -12, tau = 0.004, boson_dim = 12;
mQ1 = [0; 0.1928049391966899;
0.698712502731356;
1.3598994689006094;
2.0362294965342893;
2.6482761444031553;
3.164786931758567;
3.5905879109045467;
3.9340086185003056;
4.209022013531566;
4.4279501166106785;
4.604639179579112;
4.749590978166249;
4.861654201304321;
4.963128041625742;
5.014036022988648;
5.085586171178299;
5.133925566983976;
5.151450707226143;
]
vQ1 = [0; 2.9096781386611985;
10.533873748754877;
19.376249385981115;
27.218286001972196;
33.67099796078335;
38.14476917582421;
41.84118408443238;
44.005824732495235;
45.82748059907051;
46.92122542085746;
48.11985712473362;
49.38030241677086;
49.510725367919036;
50.33140120949116;
49.297312020870486;
50.12905662276805;
50.327785229517026;
50.4660406029108;
]
#


#T = 0.1, alpha = 1.0, N_chain = 160, maxdim = 80, cutoff = -12, tau = 0.004, boson_dim = 12;
mQ2 = [0.04950881330431919;
0.19280396896564553;
0.4157632280669486;
0.6986991417266992;
1.0200794923177707;
1.3599402450141795;
1.7023056449725407;
2.0358525793389948;
2.3541168461130715;
2.647953221183166;
2.91713024154646;
3.1633179964411293;
3.386263195790718;
3.5869280436748063;
3.768167633620232;
3.930771649265602;
4.076431819206363;
4.202329092937693;
4.315886755858702;
4.4196043505088625;
4.5122872518717125;
4.594842166535745;
4.6678011574066165;
4.736537907377009;
4.791973030277638;
4.845187341343349;
4.893794463713372;
4.944065102089293;
4.98559203808391;
5.007578521596513;
5.0437365738578475;
5.069481698611916;
5.089977912424944;
5.10762724088963;
5.123840674150535;
5.144272977843222;
5.161352869747592;
5.172047912285402;
5.181967674290515;
]
vQ2 = [0.7445004610401575;
2.9071565891006808;
6.335863222368894;
10.530655426312453;
14.608082442883937;
19.38303505415417;
23.42518887544912;
27.23287346140589;
30.745711688055014;
33.64564275635263;
36.04954625196238;
38.16346450750784;
40.0707428515116;
41.606142613212576;
43.0439379451837;
44.033764758890094;
45.09080966076878;
45.55058890859446;
46.25984736896028;
46.76202934329706;
47.342410853659175;
47.864782513433454;
48.403052398709185;
49.59631155726811;
49.46215985803795;
49.294080914545084;
49.98778970420428;
50.61210193692912;
50.810773453549054;
50.094110941007486;
50.41738562544923;
50.53948239956502;
50.30700092195174;
50.519405759565295;
50.601929832125194;
50.60791086650619;
50.98520107067508;
50.94132339252706;
50.77397626029705;
]
#


mQ1 = mQ1[1:19]
vQ1 = vQ1[1:19]
mQ2 = mQ2[1:19]
vQ2 = vQ2[1:19]


# Time step duration
tau = 0.005
nt = 250
ttotal = nt * tau  # Total time evolution
#time_steps = collect(0:10*tau:ttotal)

time_steps = collect(0:lastindex(mQ2)-1) * 10 * tau

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
xtickfont = font(12), ytickfont = font(12), ylims = (0, 6), xticks = (xticks, xtick_labels),# yticks = (yticks1, ytick_labels1), 
xguidefontsize=18,yguidefontsize=20, color = :orange,
seriesalpha = 0.9, linewidth = 6, label = L"\ \Delta = 0, T = 1.0")
plot!(time_steps, mQ2, ylims = (0, 6), color = :brown,seriesalpha = 0.9, linewidth = 6, label = L"\ \Delta = 1, T = 0.1")


plot!(twinx(),time_steps, vQ2, line = :dot, color = :brown,ytickfont = font(12),
seriesalpha = 0.9, linewidth = 6, ylims = (0, 50), label = false,yguidefontsize=20, yaxis = L"0.1\times\langle\langle Q^2 \rangle\rangle")
plot!(time_steps, vQ1, line = :dot, ylims = (0, 50), color = :orange, seriesalpha = 0.9, linewidth = 6, label = false)

plot!(widen = false)
vline!([xlims(p1)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p1)[2]], lc = :black, lw = 2, label = false)
#plot!(legend = :outertopright)

title!(L"\epsilon = 1, \alpha = 1")

#ticks_length!(tl=.04)

# Customize plot appearance, place legend inside
plot!(legend = :right, grid = false, legendfontsize = 16)

display("image/png", p1)
