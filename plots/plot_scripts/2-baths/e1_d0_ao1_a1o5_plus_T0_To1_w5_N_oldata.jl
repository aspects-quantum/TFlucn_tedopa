using DrWatson, Plots, LaTeXStrings
using Plots.PlotMeasures, QuadGK
#using PyPlot

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

# T = 0.1, alpha = 0.25, N_chain = 180, maxdim = 90, cutoff = -13, tau = 0.002, jump = 10, boson_dim = 9, omega = 5;
mQ1 = [3.868071495559938e-7;
0.3454818942048034;
0.642487154092873;
0.9770087143047327;
1.2976960323171631;
1.6007647058973404;
1.8831916632081447;
2.142894311684999;
2.3786071353609413;
2.5898631695751737;
2.7768857190339946;
2.940461224533664;
3.0818175323274297;
3.202490057360608;
3.304215016741745;
3.3888252034822437;
3.4581700346883513;
3.5140463091790717;
3.5581719569107726;
3.5920723005810853;
3.6172177222087543;
3.6349292395962873;
3.6464498997742596;
3.652867614123756;
3.655147406532976;
3.6539518283221737;
3.650028698423969;
3.643061695809405;
3.6362085123266845;
3.627946698999555;
3.618641109768604;
3.608571626811371;
3.5982237896738924;
3.5862814455722383;
3.5770455941036503;
3.567247057997769;
3.557465460700677;
3.5476844698227703;
3.5386104102963847;
3.529638841107027;
3.522662139052652;
]

vQ1 = [-1.3044521747947708e-9;
14.514909450945973;
28.178555363895008;
41.77625287161835;
54.46919218497084;
66.03978148089898;
76.33819814468154;
85.27267176483943;
92.81807355542568;
99.01252668688718;
103.93477901645649;
107.76680581806431;
110.44925920842364;
112.14001483886419;
113.12954748210976;
113.55557474385205;
113.51336178282597;
113.25731797286133;
111.97630996097665;
110.98312124498877;
109.92932144857366;
108.82879985147619;
107.3741257490419;
105.98027326202639;
104.57265344430745;
103.97186318437407;
103.22507008136648;
102.06346980288633;
100.15046512776324;
98.49994538406978;
98.23170890087563;
98.26365748789054;
97.98115204228978;
97.93742834228142;
96.94391119548402;
95.82887301452173;
96.6205061107594;
97.48293430873754;
96.6507525110492;
94.55292605743126;
93.37778083824207;
]

# T = 0.1, alpha = 1.25, N_chain = 160, maxdim = 100, cutoff = -12, tau = 0.001, jump = 20, boson_dim = 14, omega = 5;
mQ2 = [9.442235164101228e-7;
0.32618266230413073;
0.4081172432306151;
0.738322311744077;
1.0479537052931647;
1.3303380148232629;
1.5809188058793255;
1.7970334011220983;
1.9777707037127965;
2.1237572565848244;
2.236870661634194;
2.319909692402826;
2.3763445584614833;
2.4099706394280465;
2.4243623383923465;
2.423446625894705;
2.410723373677769;
2.3889704277615587;
2.3621844840642;
2.331746409105968;
2.2985039943858143;
2.264395096592755;
2.233952034528676;
2.2055633426149246;
2.180431710011721;
2.1585394974117134;
2.1378590833660205;
2.1180009391503476;
2.1092827933693106;
2.102563718833505;
]

vQ2 = [1.703886025403454e-9;
53.42178182552183;
102.30541666364284;
152.44647180782965;
199.13792810408412;
241.53801139071095;
279.0731006260495;
311.40673630935584;
338.4664039904107;
360.36627363368694;
377.4507114525304;
390.21066027297655;
398.690506061242;
404.1051175314249;
406.33526357085736;
406.9519647479768;
404.99355038630625;
401.18847709000664;
396.97389735149227;
392.25013330744423;
386.7003082139081;
382.54378833895765;
376.0775546306587;
369.5491043386846;
365.36921389265984;
361.53725331360414;
356.7001885574693;
351.16946530203853;
346.01916716531554;
340.16559426908947;
]

# T = 0, alpha = 0.25, N_chain = 180, maxdim = 90, cutoff = -12, tau = 0.002, boson_dim = 11, omega_C = 5, jump = 10;

mQ3 = [1.0060586630426573e-6;
1.7369560002639735;
3.3972551773614303;
5.074660054894658;
6.690854220829325;
8.22938982368157;
9.677111872221904;
11.024280546696238;
12.264634671218202;
13.395096418170054;
14.415483469868638;
15.328021666470603;
16.136871533647472;
16.84763051682614;
17.466892427321064;
18.001862163653193;
18.459982796127107;
18.848665212075243;
19.175069873117913;
19.446028563404425;
19.66797623563294;
19.846912600502666;
19.988397227925585;
20.097325909687285;
20.178186248459994;
20.23447172601538;
20.2697269070965;
20.286636789084017;
20.28833290913817;
20.27722525116693;
20.255969327699347;
20.226027811730827;
20.18902011940036;
]

vQ3 = [-6.390208558516089e-10;
33.546351593501065;
65.95975945520962;
97.56374463140726;
127.31577662201454;
154.76086360940909;
179.56878668015673;
201.5153226987721;
220.52488998945714;
236.62532784899446;
249.94318015163955;
260.6740370145865;
268.79527941137223;
274.9267695168528;
279.28259735764595;
282.4678523114356;
283.7442014316625;
283.64143195775387;
283.2568315839615;
282.4677164636474;
281.0628102791839;
280.24109379805924;
277.5088781219608;
274.51868140978064;
273.7021808962697;
273.4113568141816;
270.65783691986053;
267.3259836678518;
264.7736175468823;
262.1144334603907;
263.2488713534431;
267.1482573687355;
263.4668627591864;
]


#

# T = 0, alpha = 1.25, N_chain = 180, maxdim = 90, cutoff = -12, tau = 0.002, boson_dim = 14, omega_C = 5, jump = 10;
mQ4 = [5.569782460060408e-7;
1.7274362206572524;
3.2124317078215303;
4.885151963109806;
6.489014621315366;
8.00524406994917;
9.41937622643432;
10.72139565108959;
11.905625332587809;
12.970315387471485;
13.91709456704926;
14.750333738079798;
15.476432200096948;
16.10318601117009;
16.639217624475464;
17.093498716271544;
17.474956901861166;
17.7920655502673;
18.05291170720224;
18.263985000832307;
18.43491720760919;
18.569521219384495;
18.672854773482147;
18.750434484066485;
18.80617998129192;
18.84233994101655;
18.86562482465932;
18.875462886513784;
]

#T = 0, alpha = 1.25, N_chain = 190, maxdim = 90, cutoff = -11, tau = 0.001, boson_dim = 15, omega_C = 5, jump = 20;
vQ4 = [-2.9019302392483174e-9;
72.62808790851167;
141.2715964132728;
210.18455308253854;
275.3862640124725;
335.9565527685841;
391.21256746194734;
440.6968887773008;
484.2401453145255;
521.9162383307382;
553.8426949865343;
580.1903515502908;
601.643392546475;
618.7415673137098;
632.0243566147345;
642.0480934096801;
647.9628704603479;
652.3750806110804;
654.9872251948736;
657.1438793841494;
655.6878873243915;
652.8476226647033;
651.7013962073374;
650.9067385007858;
648.4822245438878;
645.911339412472;
640.6294830771567;
635.5531408893002;
]




#= 
num = 48
mQ1 = mQ1[1:num]
vQ1 = vQ1[1:num]
mQ2 = mQ2[1:num]
vQ2 = vQ2[1:num]
Sx1 = Sx1[1:num]
Sx2 = Sx2[1:num]

 =#
 ω_C = 5
 f(t, α) = α*ω_C*(t)^2 /(1+t^2)
 #= 
 g(ω, t, β) = (1-cos(ω*t)*coth(β*ω/2)*2*α*ω*exp(-ω/ω_C))
 vQ(t, α, β) = 0.5 .* [quadgk(g(ω, t[i], β), 0, 10^7) for i in 1:lastindex(t)] =#
 #= 
 g(ω, t, β, α) = (1 - cos(ω * t)) * coth(β * ω / 2) * 2 * α * ω * exp(-ω / ω_C)
 vQ(t, α, β) = 0.5 .* [quadgk(ω -> g(ω, t[i], β, α), 0, 10^7)[2] for i in 1:lastindex(t)]
  =#
 
 #= 
 function vQ(t, α, β, ω_C=5)
     g(ω, t) = (1 - cos(ω * t/ω_C)) * coth(β * ω / 2) * 2 * α * ω * exp(-ω / ω_C)
     return 0.5 .* [quadgk(ω -> g(ω, t[i]), 0, 1e5)[1] for i in eachindex(t)]
 end
  =#
 time_steps = collect(0:250) * 0.02 * ω_C
#= 
 
 mQ1_exact = f.(time_steps, 0.25)
 mQ2_exact = f.(time_steps, 1.25)
 mQ1_T_exact = f.(time_steps, 0.25)
 mQ2_T_exact = f.(time_steps, 1.25)
 
 vQ1_exact = vQ(time_steps, 0.25, 100000)
 vQ2_exact = vQ(time_steps, 1.25, 100000)
 vQ1_T_exact = vQ(time_steps, .25, 1)
 vQ2_T_exact = vQ(time_steps, 1.25, 1) =#

#time_steps4 = collect(0:lastindex(Sx1_T)-1)*10*tau*ω_C


# Set tick positions and convert tick labels to LaTeX strings automatically
xticks1 = (2:6:14) #range(2, stop = 15, length = 3)
xticks2 = (2:6:14)
xticks3 = (2:6:15) #range(5, stop = maximum(time_steps[end-1]), length = 3)
xticks3sub = (10:4:15) #range(10, stop = maximum(time_steps[end-1]), length = 2)

yticks1 = (0:3:30) #range(0, stop = 6, length = 3)
yticks2 = (0:15:650) #range(0, stop = 30, length = 3)
yticks3 = (5:20:200) #range(4, stop = 10, length = 4)
yticks3sub = range(4.5, stop = 5.5, length = 2)

#xtick_labels1 = [string(round(x, digits = 1)) for x in xticks1]
    xtick_labels1 = [" " for x in xticks1]

xtick_labels2 = [string(Int(round(x, digits = 1))) for x in xticks2]
xtick_labels3 = [string(Int(round(x, digits = 1))) for x in xticks3]
xtick_labels3sub = [string(Int(round(x, digits = 1))) for x in xticks3sub]
ytick_labels1 = [string(Int(round(y, digits = 1))) for y in yticks1]
ytick_labels2 = [string(Int(round(y, digits = 1))) for y in yticks2]
ytick_labels3 = [string(Int(round(y, digits = 1))) for y in yticks3]
ytick_labels3sub = [string((round(y, digits = 1))) for y in yticks3sub]


# Plot
gap = 2
markersize = 3.5
p1 = scatter(time_steps[1:gap:length(mQ1)], mQ1[1:gap:end], markersize=markersize, markerstrokewidth=.8, xaxis = "", 
xticks = (xticks1, xtick_labels1), yticks = (yticks1, ytick_labels1),
xtickfont = font(14), ytickfont = font(14),
xguidefontsize = 22, yguidefontsize = 17, color = :lightblue) 
scatter!(time_steps[1:gap:length(mQ2)], mQ2[1:gap:end],markersize=markersize, markerstrokewidth=.8,  color = :teal) 
scatter!(time_steps[4:gap:length(mQ3)], mQ3[4:gap:end], markersize=markersize, markerstrokewidth=.8,  color = :lightcoral)
scatter!(time_steps[4:gap:length(mQ4)], mQ4[4:gap:end], markersize=markersize, markerstrokewidth=.8, color = :brown)
plot!(yaxis = L"⟨Q⟩")
plot!(legend = false)
xlims!(0, 14)  
ylims!(0, 30)  
plot!(widen = false)
vline!([xlims(p1)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p1)[2]], lc = :black, lw = 2, label = false)
plot!(grid = false)
plot!(bottom_margin = -8mm)
annotate!(13.3, 5.2, text("(a)", 12, :black, :right))
ticks_length!(tl=.03)


gap = 2
markersize = 3.3
p2 = scatter(time_steps[1:gap:length(vQ1)], vQ1[1:gap:end], markersize=markersize, markerstrokewidth=.6,  xaxis = "", 
xticks = (xticks2, xtick_labels2), yticks = (yticks2, ytick_labels2),
xtickfont = font(14), ytickfont = font(14),
xguidefontsize = 20, yguidefontsize = 17, color = :lightblue) 
scatter!(time_steps[1:gap:length(vQ2)], vQ2[1:gap:end],markersize=markersize, markerstrokewidth=.6, color = :teal) 
scatter!(time_steps[5:gap:length(vQ3)], vQ3[5:gap:end], markersize=markersize, markerstrokewidth=.6,  color = :lightcoral)
scatter!(time_steps[5:gap:length(vQ4)], vQ4[5:gap:end], markersize=markersize, markerstrokewidth=.6,  color = :brown)
plot!(widen = false)
xlims!(0, 14)  
ylims!(0, 650)  
plot!(yaxis = L"⟨⟨Q^2⟩⟩")
plot!(xaxislabel = false)
vline!([xlims(p2)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p2)[2]], lc = :black, lw = 2, label = false)
plot!(grid = false)
plot!(legend = false)
annotate!(13.3, 26, text("(b)", 12, :black, :right))
ticks_length!(tl=.03)
plot!(xaxis = L"tω_C")

gap = 6
markersize = 3.3


p3 = scatter([], [],
    label = L"\ \ 0.0 \ \ 0.25",
    markersize = 1,     # 👈 different marker size!
    color = :lightblue, markerstrokewidth=.2
)

scatter!([], [],
    label = L"\ \ 0.0 \ \ 1.25",
    markersize = 1,     # 👈 different marker size!
    color = :teal, markerstrokewidth=.2
)

scatter!([], [],
    label = L"\ \ 1.0 \ \ 0.25",
    markersize = 1,     # 👈 different marker size!
    color = :lightcoral, markerstrokewidth=.2
)

scatter!([], [],
    label = L"\ \ 1.0 \ \ 1.25",
    markersize = 1,     # 👈 different marker size!
    color = :brown, markerstrokewidth=.2
)
#plot!(time_steps[10:end], vQ2_exact[10:end]./mQ2_exact[10:end], color = :teal, seriesalpha = 0.6, linewidth = 2, label = false)
#plot!(time_steps[10:end], vQ2_T_exact[10:end]./mQ2_T_exact[10:end], color = :brown, seriesalpha = 0.6, linewidth = 2, label = false)

scatter!(time_steps[1:gap:length(vQ1)], vQ1[1:gap:end] ./ mQ1[1:gap:end], color = :lightblue, xticks = (xticks3, xtick_labels3), yticks = (yticks3, ytick_labels3), xtickfont = font(14), ytickfont = font(14),
	xguidefontsize = 20, yguidefontsize = 20, seriesalpha = 1, xaxis = L"tω_C", markersize=markersize, label = "") #, label = L"\  (1, 1.5)")
scatter!(time_steps[2:gap:length(vQ2)], vQ2[2:gap:end] ./ mQ2[2:gap:end], markersize=markersize, color = :teal, label = "")
scatter!(time_steps[3:gap:length(vQ3)], vQ3[3:gap:end] ./ mQ3[3:gap:end], markersize=markersize, markerstrokewidth=.6, color = :lightcoral, label = "")
scatter!(time_steps[4:gap:length(vQ4)], vQ4[4:gap:end] ./ mQ4[4:gap:end],markersize=markersize, markerstrokewidth=.6, color = :brown, label = "")
plot!(widen = false)
ylims!(4, 200)  
xlims!(0, 14)  
vline!([xlims(p3)[1]], lc = :black, lw = 2, label = false)
hline!([ylims(p3)[2]], lc = :black, lw = 2, label = false)
plot!(grid = false, ymirror = true)
#plot!(legend = false)
annotate!(16.5, 7, text(L"F", 18, :black, :right))
ticks_length!(tl=.02)
annotate!(13.4, 9.7, text("(c)", 12, :black, :right))
plot!(legendtitle=L"\ \ \ \ T \ \ \ \ \ \alpha", legendtitlefontsize=12,legend = (.6,.83),legendfontsize = 12)
plot!(foreground_color_legend = RGBA(0, 0, 0, 0.1), background_color_legend =  RGBA(0, 0, 0, 0))
#= plot!(legend = true)
#plot!(legend = :topright)
plot!(legend = (.55,.8), grid = false, legendfontsize = 10, legend_marker_size = 0)
plot!(legendtitle=L"\ \ \ \ T \ \ \ \ \ \alpha", legendtitlefontsize=11)
plot!(foreground_color_legend = RGBA(0, 0, 0, 0.1), background_color_legend =  RGBA(0, 0, 0, 0)) =#


#
#= 
plot!(p3, inset=bbox(0.35,0.01,0.42, 0.6), subplot=2)
gap = 5
markersize = 4
scatter!(p3[2], time_steps[50:gap:length(vQ1)], vQ1[50:gap:end] ./ mQ1[50:gap:end], color = :lightblue, xticks = (xticks3sub, xtick_labels3sub), yticks = (yticks3sub, ytick_labels3sub), xtickfont = font(12), ytickfont = font(12),
	xguidefontsize = 10, yguidefontsize = 15, seriesalpha = 1, xaxis = "", markersize=markersize, markerstrokewidth=.6,) #, label = L"\  (1, 1.5)")
scatter!(p3[2], time_steps[50:gap:length(vQ2)], vQ2[50:gap:end] ./ mQ2[50:gap:end], markersize=markersize, markerstrokewidth=.6, color = :teal) #, label = L"\ (0.1, 0.75)")
scatter!(p3[2], time_steps[50:gap:length(vQ1_T)], vQ1_T[50:gap:end] ./ mQ1_T[50:gap:end], markersize=markersize, markerstrokewidth=.6, color = :lightcoral)
scatter!(p3[2], time_steps[50:gap:length(vQ2_T)], vQ2_T[50:gap:end] ./ mQ2_T[50:gap:end], markersize=markersize, markerstrokewidth=.6, color = :brown)
ylims!(p3[2], 4.4, 5.9)  
xlims!(p3[2], 7, 14)  
plot!(p3[2], grid = false)
plot!(p3[2], legend = false, aspect_ratio = 5., framestyle = :box)	 =#



custom_layout = @layout [[a{0.5h}; b{1.35w}] c{0.6w}]
p = plot(p1, p2, p3, layout = custom_layout, size = (600, 400), left_margin = 5mm, right_margin = 9mm, dpi=600)



#p=plot(p1, p2, p3, p4, layout=(2,2))
display("image/png", p)
savefig(p, "plot_4.png")