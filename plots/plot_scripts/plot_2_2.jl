

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


T = 0.1
α = 0.5

#T = 1, alpha = 1.0, N_chain = 160, maxdim = 80, cutoff = -12, tau = 0.004, boson_dim = 12;

mQ1 = [0; 0.1928065160742508;
0.6980226449507666;
1.3585007147905535;
2.0333771468369815;
2.645421586208445;
3.1595761150466353;
3.579991543886029;
3.919323288250631;
4.183997916991227;
4.404580409216716;
4.566594295654851;
4.694059588076129;
4.804921117964087;
4.8894182710300464;
4.962157240937958;
]

vQ1 = [0; 2.91672005054159;
10.586712048634627;
19.487261862492943;
27.262565137914347;
33.36032163498322;
38.09661063471608;
41.54902428481861;
43.5484252216949;
45.00017758500818;
47.31101201531439;
48.03272038409743;
49.26478107001513;
50.07193325344016;
50.776924265175296;
50.80751240161593;
]


#


#T = 0.1, alpha = 1.0, N_chain = 160, maxdim = 80, cutoff = -12, tau = 0.004, boson_dim = 12;
mQ2 = [0; 0.1928049391966899;
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
vQ2 = [0; 2.9096781386611985;
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

#= 
mQ1 = mQ1[1:16]
vQ1 = vQ1[1:16]
mQ2 = mQ2[1:16]
vQ2 = vQ2[1:16]

 =#
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
xtickfont = font(12), ytickfont = font(12),
xguidefontsize=18,yguidefontsize=20, color = :orange,
seriesalpha = 0.9, linewidth = 6, label = L"\  T = 1.0")
plot!(time_steps, mQ2, color = :brown,seriesalpha = 0.9, linewidth = 6, label = L"\ T = 0.1")


plot!(time_steps, vQ1, line = :dot, color = :orange, seriesalpha = 0.9, linewidth = 6, label = false)
plot!(twinx(),time_steps, vQ2, line = :dot, color = :brown,ytickfont = font(12),
seriesalpha = 0.9, linewidth = 6, label = false,yguidefontsize=20, yaxis = L"0.1\times\langle\langle Q^2 \rangle\rangle")

plot!(widen = false)
vline!([xlims(p1)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p1)[2]], lc = :black, lw = 2, label = false)
#plot!(legend = :outertopright)

title!(L"\epsilon = 1, \Delta = 0, \alpha = 1")

#ticks_length!(tl=.04)

# Customize plot appearance, place legend inside
plot!(legend = :bottomright, grid = false, legendfontsize = 14)

display("image/png", p1)
