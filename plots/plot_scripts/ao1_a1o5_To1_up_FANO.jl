

using DrWatson, Plots, LaTeXStrings
using Plots.PlotMeasures

#GR()  # Use GR backend, but you can switch to PyPlot if needed
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


#T = 0.1, alpha = 0.1, N_chain = 160, maxdim = 80, cutoff = -11, tau = 0.002, jump = 10, boson_dim = 12, omega = 10;

mQ1 = [-8.962948512291943e-10;
0.03852094318899247;
0.1389126362277033;
0.26880223012714977;
0.40026071619466347;
0.5184501303848668;
0.6188269630276709;
0.7020696122507838;
0.7706915070400471;
0.82745803893244;
0.8747995680720427;
0.9146873108347314;
0.9486625602814273;
0.9779145497792251;
1.0033516425366475;
1.025675711326448;
1.0454297297455886;
1.0630596194539041;
1.0789152949845686;
1.093276559038289;
1.106374437431584;
1.118398591452681;
1.1295118931972472;
1.1398436212292669;
1.1495040667644405;
1.1585788723151529;
1.1671387918958402;
1.175236641388366;
1.1829223520209484;
1.1902384610731802;
1.1972104811954356;
1.2038786876185654;
1.2102630699922012;
1.2163876560865037;
1.2222769512905374;
1.227935156139231;
1.2333805654295582;
1.2386387637713525;
1.2437139882987813;

]
vQ1 = [6.851946888353368e-19;
0.0014849583703256267;
0.019310172010917594;
0.07230222904437908;
0.1603085266662633;
0.2689498894276146;
0.3831631127445231;
0.49316904251966304;
0.5942792622051858;
0.6850439310066662;
0.7656728841538328;
0.8370899840457945;
0.9004361509837286;
0.9568312906584684;
1.0072713559980089;
1.0526095292481772;
1.0935704989044648;
1.1307915248474403;
1.164804571096725;
1.1960551779420294;
1.224921379169047;
1.2517321403545865;
1.2767742107501439;
1.3002814087383812;
1.3224574315318471;
1.3434628606733072;
1.3634271164066503;
1.382455491625455;
1.4006424171385305;
1.4180662218466853;
1.4347844960269691;
1.4508660354788523;
1.4663541192421548;
1.4812972131864597;
1.4957334720786926;
1.5096803355563704;
1.5231727845545988;
1.5362589057758735;
1.548947767524782;
]
#


# T = 0.1, alpha = 1.5, N_chain = 160, maxdim = 80, cutoff = -11, tau = 0.002, jump = 10, boson_dim = 12, omega = 10;
mQ2 = [-8.038587013413791e-10;
0.5768138100246434;
2.08038459691493;
4.019933084382826;
5.950482710766394;
7.69228073902944;
9.122880406336707;
10.266856492548536;
11.168020258580004;
11.879985300526222;
12.435347305343505;
12.853196002348938;
13.245908043671403;
13.502538515043607;
13.762987305195672;
14.001766840471287;
14.161516198111965;
14.266860014827587;
]
vQ2 = [5.620333592929843e-19;
16.39912112170579;
63.80781012401746;
118.97008960745033;
172.58264029189925;
215.7687796522311;
242.38908990891977;
266.6726682108764;
288.4221766936284;
297.2364688008148;
307.0743469055973;
315.08546113634924;
322.10545405142386;
322.56058661431933;
326.654387552223;
325.31557980934303;
329.46723021045597;
327.6456250203603;
]
#


mQ1 = mQ1[1:18]
vQ1 = vQ1[1:18]
mQ2 = mQ2[1:18]
vQ2 = vQ2[1:18]


# Time step duration
tau = 0.002
nt = 250
omega_C = 10
ttotal = nt * tau  # Total time evolution
#time_steps = collect(0:10*tau:ttotal)

time_steps = collect(0:lastindex(mQ2)-1) * 10 * tau *omega_C

# Set tick positions and convert tick labels to LaTeX strings automatically
xticks = range(10 * tau *omega_C, stop = maximum(time_steps), length = 5)
yticks1 = range(0, stop = 6, length = 5)
yticks2 = range(0, stop = 60, length = 5)

xtick_labels = [string(round(x, digits = 1)) for x in xticks]
ytick_labels1 = [string(round(y, digits = 1)) for y in yticks1]
ytick_labels2 = [string(round(y, digits = 1)) for y in yticks2]


p1 = plot(time_steps, vQ1./mQ1, xaxis = L"tω_C", yaxis = L"F", #yaxis = L"\frac{\langle\langle Q^2 \rangle\rangle}{\langle Q \rangle}", 
xtickfont = font(20), ytickfont = font(20), xticks = (xticks, xtick_labels),
xguidefontsize=28,yguidefontsize=28, color = :lightcoral,
seriesalpha = 1, linewidth = 6, label = "")

# Add a dummy dataset with markers for the legend
#= scatter!(time_steps, vQ1 ./ mQ1,
      color = :cornflowerblue, marker = :circle, markersize=4, markerstrokewidth=0,# Marker-only for legend
      label = L"\ \alpha = 0.1")

# Add a dummy dataset with markers for the legend
scatter!(time_steps, vQ2 ./ mQ2, seriesalpha=0.5,
      color = :violetred4, marker = :circle, markersize=4, markerstrokewidth=0, # Marker-only for legend
      label = L"\ \alpha = 1.5") =#
plot!(time_steps, vQ2./mQ2, color = :teal,seriesalpha = .8, linewidth = 6, label = "")




plot!(widen = false)
vline!([xlims(p1)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p1)[2]], lc = :black, lw = 2, label = false)
#plot!(legend = :outertopright)
#annotate!(-0.1, 3, Plots.text(L"Force (F)", 15, :black, rotation=90))
#title!(L"\epsilon = 1, \Delta = 0, T = 0.01")
plot!(size=(500, 350))

ticks_length!(tl=.02)


# Customize plot appearance, place legend inside
plot!(legend = :right, grid = false, legendfontsize = 20, left_margin = [2mm 0mm], right_margin = [10mm 10mm])   #,framestyle = :box), legendposition = (1.6, 0.7)

display("image/png", p1)
