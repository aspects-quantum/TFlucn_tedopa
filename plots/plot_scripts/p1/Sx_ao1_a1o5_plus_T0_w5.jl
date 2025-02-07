

using DrWatson, Plots, LaTeXStrings
using Plots.PlotMeasures

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


#T = 0, alpha = 1.5, N_chain = 160, maxdim = 80, cutoff = -12, tau = 0.005, jump = 5, omega_C = 5, boson_dim = 13;
Sx2 = [0.5000000000000079;
0.4999937539153088;
0.4997788975180781;
0.4992608400104756;
0.4984682694772785;
0.4974238988186841;
0.4961715698239579;
0.49474481233133827;
0.4931870233412073;
0.49153918159946874;
0.48984072023888736;
0.4881272217767965;
0.48642747101466566;
0.4847715979701658;
0.4831820321147927;
0.4816713648177493;
0.48025288020973655;
0.47893338744139;
0.47771379867712604;
0.476594339459994;
0.47557395791521373;
0.47464589446634375;
0.47380551114060365;
0.473047023539649;
0.4723633889130105;
0.47174972805520154;
0.47120129084630585;
0.4707091044035708;
0.4702690286716751;
0.4698745028980611;
0.4695221171607096;
0.46920673923851286;
0.4689245135551995;
0.46867209546444905;
0.46844539748823727;
0.4682415266568821;
0.46805912406845707;
0.467895168019274;
]



#T = 0, alpha = 0.1, N_chain = 160, maxdim = 80, cutoff = -12, tau = 0.005, jump = 5, omega_C = 5, boson_dim = 13;
Sx1 = [0.5000000000000079;
0.499993750273277;
0.4997753516538661;
0.499247653096056;
0.49841215415803913;
0.4972661512097529;
0.49581932107717225;
0.49407537212002944;
0.49203325766001244;
0.4897056074327806;
0.48709013164392323;
0.48419879381322284;
0.48103043354971153;
0.4775961598862138;
0.4738982850171625;
0.4699435127729113;
0.4657373069311448;
0.4612854857326785;
0.4565938542497714;
0.451668895513928;
0.4465157011281731;
0.4411405980976286;
0.435546069028752;
0.4297388637015962;
0.42372427396143;
0.4175082592797544;
0.411098999754206;
0.4044972362533452;
0.3977100242902663;
0.3907437503133952;
0.38360461563657783;
0.37629809354583876;
0.36882743840349347;
0.3611992431494381;
0.35341994792935555;
0.3454965730669159;
0.3374304621156686;
0.3292300902830357;
0.3209005502341082;
0.3124482193220841;
0.30387770964382455;
0.2951932040502059;
0.2864019096403889;
0.2775096705818049;
0.26852126072708177;
0.25944147751261026;
0.25027736848380194;
0.2410339688932861;
0.2317162313187679;
0.2223304280611705;
0.2128813619366386;
0.20337505751022905;
0.1938164600388608;
0.18421116840602506;
0.17456473650111734;
0.16488242755275376;
0.1551694681631056;
0.14543098532796755;
0.13567265602036144;
0.12589976327143662;
0.11611714776324002;
0.10632984014850906;
0.09654369638667143;
0.0867633433789585;
0.07699374900441934;
0.06724042572013485;
0.05750775620196002;
0.047800916502853015;
0.038124092997061584;
0.028482476069063272;
0.0188812846675323;
0.009325600778337042;
-0.00018071170622296085;
-0.009633175598690875;
-0.01902662945980023;
-0.02835671315382775;
-0.03761930117885838;
-0.04681023818251692;
-0.0559255593238911;
-0.064960239263983;
-0.07391102886450211;]
#


Sx1 = Sx1[1:38]
Sx2 = Sx2[1:38]

# Time step duration
tau = 0.005
omega_C = 5
#time_steps = collect(0:10*tau:ttotal)

time_steps = collect(0:lastindex(Sx2)-1) * 5 * tau*omega_C

# Set tick positions and convert tick labels to LaTeX strings automatically
xticks = range(0, stop = maximum(time_steps), length = 5)
yticks1 = range(0, stop = 6, length = 5)
yticks2 = range(0, stop = 60, length = 5)
#= 
xtick_labels = [L"\$" * string(round(x, digits = 2)) * "\$" for x in xticks]
ytick_labels1 = [L"\$" * string(round(y, digits = 1)) * "\$" for y in yticks1]
ytick_labels2 = [L"\$" * string(round(y, digits = 1)) * "\$" for y in yticks2] =#

xtick_labels = [string(round(x, digits = 1)) for x in xticks]
ytick_labels1 = [string(round(y, digits = 1)) for y in yticks1]
ytick_labels2 = [string(round(y, digits = 1)) for y in yticks2]


p1 = plot(time_steps, Sx1, xaxis = L"tω_C", yaxis = L"S_x", #yaxis = L"\frac{\langle\langle Q^2 \rangle\rangle}{\langle Q \rangle}", 
xtickfont = font(20), ytickfont = font(20), xticks = (xticks, xtick_labels),
xguidefontsize=28,yguidefontsize=28, color = :teal,
seriesalpha = 1, linewidth = 9, label = "")

# Add a dummy dataset with markers for the legend
#= scatter!(time_steps, Sz1,
      color = :cornflowerblue, marker = :circle, markersize=5, markerstrokewidth=0,# Marker-only for legend
      label = L"\ \alpha = 0.1")

# Add a dummy dataset with markers for the legend
scatter!(time_steps, Sz2, seriesalpha=0.5,
      color = :violetred4, marker = :circle, markersize=5, markerstrokewidth=0, # Marker-only for legend
      label = L"\ \alpha = 1.5") =#
plot!(time_steps, Sx2, color = :lightcoral,seriesalpha = 1, linewidth = 9, label = "")


plot!(widen = false)
vline!([xlims(p1)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p1)[2]], lc = :black, lw = 2, label = false)
#plot!(legend = :outertopright)
plot!(size=(500, 400))

#title!(string(L"\epsilon = 1, \Delta = 0, ω_C = 10, |1\rangle"))

ticks_length!(tl=.02)

# Customize plot appearance, place legend inside
plot!(legend = false, grid = false, legendfontsize = 20, left_margin = [2mm 0mm], right_margin = [10mm 10mm])   #,framestyle = :box), legendposition = (1.6, 0.7)

display("image/png", p1)
