

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

mQ1 = [0; 0.03841005356687807;
0.14686462157712626;
0.30777945537135754;
0.49920650801011807;
0.7009170575709426;
0.8983545683864623;
1.0821859092827897;
1.2477885612668866;
1.3943517961110423;
1.5217804696887434;
1.6320489890404337;
1.7279422624575644;
1.8113179094567964;
1.8834678734040056;
1.946116324679857;
2.0015176188175676;
2.0498557385334473;
2.090750532732996;
2.127420293154871;
2.160358825977795;
2.186483675287454;
2.21098687506732;
2.2762378934288687;
2.2535328309765976;
2.2672230264452695;
2.28495873195419;
2.2992600911964027;
2.2931103096076373;
2.3342912838752303;
2.3290024273487364;
2.336808955143919;
2.3531324422921123;
2.3569476931855213;
2.366125210853806;
2.3700841109603723;
2.3742653544169547;
2.371065216711606;
2.3642221675344333;
2.3709170531584878;
2.3812254694755985;
2.3876738337008847;
2.3922736205839366;
2.398436200233744;
2.4067310431681896;
2.4180248729081715;
2.419402723159677;
2.4258293657630055;
2.426682500919342;
2.4293948762036606;
2.431151316686351;
2.4328159655395716;
2.4365468623749624;
2.4359327356879463;
2.435437536321621;
]
vQ1 = [0; 0.5671251842264959;
2.1363291792616965;
4.328919679655172;
6.721322970823623;
9.028028433283927;
11.055449350670422;
12.718389424467794;
14.049801903571206;
15.114716331868845;
15.918339161417418;
16.539400186784345;
17.034897699279036;
17.373538014869027;
17.640876267304474;
17.806620049329847;
17.984428451610817;
18.15050281950623;
18.25123194644664;
18.340150178463283;
18.36050672688766;
18.36257639885574;
18.384409151508176;
19.128192992691122;
18.349264350611413;
18.315269321527673;
18.31873764988013;
18.11669133751424;
17.754030496723807;
18.72306818682729;
17.994289489420833;
17.93413520004542;
18.021346816361355;
17.679668500324702;
17.90070298351313;
17.578326595205773;
17.43541670406494;
17.12805241084381;
17.10657768515711;
17.208635063581326;
17.659918546744844;
17.711573206627644;
17.69443750294204;
17.734341230360243;
17.53970164900784;
17.447165988097545;
17.129577715199275;
17.04630758539691;
17.122633465905224;
16.799889067221603;
16.80118322935168;
16.19297228784562;
16.17705395178781;
16.02300677307364;
16.68303583040874;
]
#



mQ1 = mQ1[1:25]
vQ1 = vQ1[1:25]

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
seriesalpha = 0.9, linewidth = 6, label = false, yguidefontsize=20, yaxis = L"0.1\times\langle\langle Q^2 \rangle\rangle")

plot!(widen = false)
vline!([xlims(p1)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p1)[2]], lc = :black, lw = 2, label = false)
#plot!(legend = :outertopright)

title!(L"\epsilon = 1, \Delta = 0, T = 0.01")

#ticks_length!(tl=.04)

# Customize plot appearance, place legend inside
plot!(legend = :right, grid = false, legendfontsize = 16)

display("image/png", p1)
