using DrWatson, Plots, LaTeXStrings
using Plots.PlotMeasures

gr()  # Use GR backend, but you can switch to PyPlot if needed
function ticks_length!(;tl=0.02)
    p = Plots.current()
    xticks, yticks = Plots.xticks(p)[1][1], Plots.yticks(p)[1][1]
    xl, yl = Plots.xlims(p), Plots.ylims(p)
    x1, y1 = zero(yticks) .+ xl[1], zero(xticks) .+ yl[1]
    sz = p.attr[:size]
    r = sz[1]/sz[2]
    dx, dy = tl*(xl[2] - xl[1]), tl*r*(yl[2] - yl[1])
    plot!([xticks xticks]', [y1 y1 .+ dy]', c=:black, labels=false)
    plot!([x1 x1 .+ dx]', [yticks yticks]', c=:black, labels=false, xlims=xl, ylims=yl)
    return Plots.current()
end
# Define data


T = 0.1
α = 0.5

# T = 1, alpha = 1.5, N_chain = 160, maxdim = 80, cutoff = -12, tau = 0.002, boson_dim = 12
#quantum_2
mQ1 = [-6.701200420323295e-10, 0.07417934231980416, 0.2882139341359709, 0.6187343643744782, 1.0337904723393054,
 1.499001940100597, 1.9839486117071081, 2.464470262387435, 2.924560812364745, 3.3538791607561578, 3.74729226071818,
  4.102875609971584, 4.422491071027604, 4.707693894350767, 4.961100636581624, 5.186737613320515, 5.3871970409006975, 
  5.564906715449066, 5.723495005840011, 5.864840574262888, 5.9906413862628565, 6.1036321345983255, 6.204957053947846, 
  6.29569604252791, 6.377948184011834, 6.452184135432528, 6.5190190599333055, 6.580048027413805, 6.635533810603567,
   6.685676426356201, 6.7318140548728636, 6.77411381766044, 6.812498267313439, 6.848117825671151, 6.88093978075531,
    6.9107604748836895, 6.93853537143811, 6.9642092657040715, 6.98776343542756, 7.009870155329719, 7.0303359149643985,
     7.049098297734362, 7.066897967572002, 7.08331545904547, 7.098482775958441, 7.11295626078314, 7.126454471034044, 
     7.138851236907995, 7.15071187230823, 7.161749227072473, 7.172028786535627, 7.181944775234197, 7.191233090533547, 
     7.199776791021262, 7.2080669247067775, 7.215848897058456, 7.223132800808743, 7.230169073441292, 7.236804308121091,
      7.242846943718516, 7.248867170026412, 7.254572039383291, 7.259943634609954, 7.264984199882683, 7.26986836445171,
       7.274428972943431, 7.278681000167391, 7.28301297074741, 7.291729582626649, 7.297365009900281, 7.303053582226106, 
       7.306606735773159, 7.309111212640251, 7.31388080831138, 7.3161180568689055, 7.318227739887866, 7.31936875657472,
        7.325117073238866, 7.3184557178448815, 7.320601160862432, 7.321610160294002, 7.327607170755406, 7.330122247268546, 
        7.3335137441795, 7.3317751651822105, 7.33314519383405, 7.338334026057966, 7.344108028666778, 7.384039349967663, 
        7.388614725897252, 7.386851537813287, 7.388633682883379, 7.388449079103263, 7.3875091022737065, 7.386567198496311,
         7.386301371118093, 7.384795920989468, 7.382384214481678, 7.380042550635595, 7.375678546045286, 7.372313401327544, 
         7.3695630332262425, 7.367425808720563, 7.365846337225915, 7.362536285152483, 7.362912983075788, 7.363318443207944,
          7.364362102900519, 7.365954996798948, 7.368510610756241, 7.371707625481829, 7.374559165212885, 7.377732531609419, 
          7.381548011917173, 7.384085239381633, 7.3860059652863095, 7.386158450162147, 7.385676020117676, 7.384558108283324, 
          7.382884882097535, 7.380598358891979, 7.378268772528979, 7.37577922129756, 7.3732332178973365, 7.373496842200399,
           7.3709323865098, 7.369980337734967, 7.368997043381827, 7.368430181320545, 7.367870727736731, 7.366992759812648,
            7.366692700127762, 7.366286492957367, 7.365749735192851, 7.365489205655934, 7.365572710566176, 7.366431118136779,
             7.364656360109717, 7.364733312439598, 7.36544017244189, 7.366741991315852, 7.36767822194565, 7.367894464677202, 
             7.367373225671076, 7.3668546894078055,
 7.366340652072987, 7.3651555827183515, 7.364715795034263, 7.363800605732896, 7.363197914567956, 7.363037212051773]

vQ1 = [3.9319661971407324e-19, 1.0885603712322123, 4.379556891289507, 9.431913888687689, 15.961497856445249, 
22.155121138188075, 28.66600310002076, 34.963431572946774, 40.846064523008586, 46.08762847701915, 50.666521957657565, 
54.51111697586523, 57.80903321319035, 61.62382743362458, 64.06642309931664, 66.12815495189663, 69.12233256696024, 
70.35104665860925, 72.67141226853498, 73.80946541058563, 75.60848818583024, 76.99358749868557, 77.73539460465953, 
78.98588879436443, 79.53853688468958, 80.44154425364687, 81.39125623801552, 82.28789389459554, 82.59580248122752,
 83.32908918780413, 83.89194160458158, 84.0805574431279, 84.67754877365773, 85.12966863737883, 85.35533295837182,
  85.96334318845267, 86.42084351311358, 86.76455851860709, 86.77134956603435, 87.17463129194262, 87.46350327537, 
  87.76390171799125, 87.68407091353139, 88.12471389102795, 88.36588942112941, 88.63137699311453, 88.59298027382476, 
  88.77740234985738, 89.00676580373297, 89.21045614693541, 89.37078140722973, 89.56410060270943, 89.49635762269932,
   89.63285742696017, 89.79318001087387, 89.92502334058562, 89.81225995003922, 89.97714668119262, 90.1158204599765,
    90.21602012692375, 90.32938573238464, 90.44511921752957, 90.36167584077083, 90.44074808674074, 90.52901377616556, 
    90.63655092750331, 90.69117995344536, 90.79289999676669, 90.92973751789386, 91.04204466801396, 91.03952415340294, 
    91.12424961771762, 91.1595162405358, 91.24725601932178, 91.23068260983104, 91.18364752627112, 91.08462512177574,
     91.01154493242704, 90.89442944260456, 91.00466866164186, 91.00300222298357, 91.29510255506773, 91.38317048371242, 
     91.62240924835397, 91.24264864600286, 91.72943810524255, 91.59619761174147, 91.81299025837646, 92.96648712416297,
      93.10475367220039, 93.15455972101952, 93.16570985589156, 93.20913619968988, 93.17232018275936, 93.15113377083097,
       93.1923281704149, 93.15883371052762, 93.19470660872221, 93.19722169849173, 93.13465511475925, 93.08398021341432, 
       93.07701021493403, 93.09356779600378, 93.08857772281435, 93.01544143715866, 93.08445514652341, 93.09594311296838,
        93.08849005947202, 93.15306264046598, 93.16876233024266, 93.21347955367622, 93.28129737778863, 93.30768503259436, 
        93.41193784151271, 93.4361923388638, 93.48826939237343, 93.50530367460347, 93.47886072884224, 93.51045200227398, 
        93.53306399501543, 93.49189016510276, 93.48795131105747, 93.41200315618313, 93.3263505115019, 93.31187011100738, 
        93.05891995163928, 92.95331745402767, 92.87830507276489, 92.82382792202817, 92.81905876875325, 92.77904233225206,
         92.8889643834166, 92.87315978478502, 92.81542939348802, 92.81446418875574, 92.76244216473674, 92.78209982063396,
          92.74838169598151, 92.75499379407253, 92.80656497679902, 92.91707334813373, 92.97336714351474, 92.97678826964712, 
          92.96224174266023, 92.91354653464684, 92.89433280243962, 92.80364147220608, 92.80193983179092, 92.76847368467489, 92.7345597311112, 92.76565216821976]


#


#T = 0.1, alpha = 0.1, N_chain = 150, maxdim = 70, cutoff = -11, tau = 0.005, boson_dim = 11;
#quantum_1
mQ2 =[-6.271436287263679e-10, 0.03708242235245901, 0.14407483608562807, 0.30930251631407624, 0.5166494011786126, 0.7490853043755736, 
0.9913929590492002, 1.2316970500724735, 1.4616800176780838, 1.6762955406513824, 1.8725915630857382, 2.0506790680549027, 2.210531356738695,
 2.3532288054077948, 2.4797864746629674, 2.592644024349136, 2.692921420178068, 2.782054885406042, 2.860947289414441, 2.9316760719078627,
  2.9948485025292184, 3.0514028410175182, 3.1017946961564555, 3.147473392744561, 3.188629443291526, 3.225799808510108, 3.2591621302348686,
   3.289783743721673, 3.3175964456208495, 3.342587076555704, 3.3658120626788506, 3.3870599308748175, 3.40649721815626, 3.42409378793723,
    3.4405953169188837, 3.4557522745980123, 3.469445750887407, 3.482493200771575, 3.4945180129173803, 3.5056422275230066, 3.515734198245793, 
3.5254661540938748, 3.5344547971224207, 3.542615356679401, 3.5504804826092746, 3.5578601063977158, 3.5647033469962413, 3.570943903931835, 3.5770326539764614, 3.58278175040808, 3.588126775179008]
vQ2 = [9.027387866221967e-20, 0.0008748509857108887, 0.01044749933013392, 0.03544259537029042, 3.662306679865938, 0.11308581382195892,
 3.075292400767091, 0.19145880087313982, 0.2258283708347382, 0.2567766074332358, 0.2845819751121633, 0.3095880638300052, 0.33211346048936125,
  0.3524662770381898, 0.3709075292516892, 0.38766044768113206, 3.567276623868815, 3.208262037197868, 0.4297001780454063, 2.7940899295612907,
   0.4523318028403185, 0.4623530910821785, 0.47162604509551537, 0.4802106452997509, 0.48817005831523436, 0.4955621377388377, 
   0.5024295924697422, 1.5738362263815906, 0.5147921635504276, 0.5203744355252127, 0.5448257036437468, 1.80681960122693, 2.239764950782718,
    1.5558629124031906, 1.7971520818463311, 2.058416690234614, 1.726230960089934, 2.558346063379364, 1.3498843227194033, 2.463174361828707, 
	1.931178209433431, 1.4464286056152529, 2.5770283547146837, 2.0692922105923204, 
0.56784714419539, 0.6437146760571453, 1.5044159081113102, 1.4516915850650471, 2.108710463141455, 2.180374332986758, 2.1731651103561167]
#

#T = 1, alpha = 1.0, N_chain = 150, maxdim = 80, cutoff = -11, tau = 0.005, boson_dim = 11;
#classical_1
mQ3 = [-4.170556387513764e-10, 0.037097508626511844, 0.1440050554060342, 0.3093191172035412, 0.5166851160161364, 0.7491641065774199, 0.9914895290418343, 
1.231768102270968, 1.461738260895286, 1.676364913641656, 1.8726876705584363, 2.0507655746563387, 2.210591560176421, 2.35324125075261, 2.4797684375179787,
 2.5925742290350633, 2.6927717905420745, 2.781851823300509, 2.8607206940079783, 2.9313935329742735, 2.9945306028930947, 3.05100704073905, 3.101337440992458,
  3.146925645785261, 3.188002883444999, 3.2246884366484596, 3.2583211238656156, 3.2888915881875755, 3.3166418228104853, 3.3415588660250957, 3.364724375527057,
   3.385838627566864, 3.4048440462110894, 3.4225809193353007, 3.4390014387798042, 3.4540450088949277, 3.4676729059625733, 3.4805455942203976, 3.4924620870125467, 
   3.5031601937889167, 3.5133880278036664, 3.5228695293995873, 3.5317487874667095, 
3.539773095433019, 3.5475140173150126, 3.5546756851887045, 3.5613365807609307, 3.567384833358774, 3.5733420047953848, 3.5788664357495676, 3.5840409907186332]

vQ3 = [1.6882132800011753e-19, 0.5673129007786027, 2.1768855568946486, 4.6772408911579735, 7.703140462232257, 11.176104827179294,
 13.770328786851923, 16.330123631379113, 18.64609537591734, 20.677434138315547, 22.345022110401068, 23.762294411028606, 24.594931055203375, 
 25.501141767158053, 26.173161231388313, 26.719090301302913, 27.153558331710027, 27.476244206809135, 27.695479072374162, 28.34760107282132,
  28.51942309833239, 28.630144455912614, 28.693132476753036, 29.29347540788001, 29.336048697477466, 29.183763523265288, 29.632690518648033, 
  29.645533761030542, 29.630320361916397, 29.951701597308933, 29.92451101087088, 30.119945237379575, 30.057346305101493, 30.291777954269175,
   30.254476118628258, 30.39015243024729, 30.31692978416687, 30.50051672165845, 30.44891896265014, 30.593293546941332, 30.483451670199855, 30.63776502237652, 
30.57600527241002, 30.703591687253574, 30.599197022747035, 30.79689912734348, 30.663798096696986, 30.77575805564748, 30.727767371892963, 30.85401664596877, 30.72275054265546]


num=50
mQ1 = mQ1[1:num]
vQ1 = vQ1[1:num]
mQ2 = mQ2[1:num]
vQ2 = vQ2[1:num]
mQ3 = mQ3[1:num]
vQ3 = vQ3[1:num]

ω_C = 5
f(t, α) = α*ω_C*(t)^2 /(1+t^2)


# Time step duration
tau = 0.002  
nt = 500
ttotal = nt * tau  # Total time evolution
#time_steps = collect(0:10*tau:ttotal)

time_steps = collect(0:lastindex(mQ3)-1)*10*tau*ω_C

mQ1_exact = f.(time_steps, 1.5)
mQ2_exact = f.(time_steps, 0.75)
mQ3_exact = f.(time_steps, 0.75)


# Set tick positions and convert tick labels to LaTeX strings automatically
xticks = (0 : 1 : maximum(time_steps[end-1]))[2:end] #range(0.2, stop = maximum(time_steps[end-1]), length = 6)
yticks1 = (0 : 7 : 7) 
yticks2 = (0 : 85 : 85) 

xtick_labels = [string(Int(round(x, digits = 1))) for x in xticks]
ytick_labels1 = [string(round(y, digits = 1)) for y in yticks1]
ytick_labels2 = [string(Int(round(y, digits = 1))) for y in yticks2]

p1 = plot(time_steps, mQ1_exact, color = :lightgray, seriesalpha = 1, linewidth = 5, label = false)
plot!(time_steps, mQ2_exact, color = :lightgray, seriesalpha = 1, linewidth = 5, label = false)
plot!(time_steps, mQ3_exact, color = :lightgray, seriesalpha = 1, linewidth = 5, label = false)
scatter!(time_steps[1:2:end], mQ1[1:2:end], markersize=10, markerstrokewidth=1, yaxis = L"\langle Q \rangle",  #xaxis = L"tω_C", 
xticks = (xticks, xtick_labels), yticks = (yticks1, ytick_labels1),
xtickfont = font(24), ytickfont = font(24), 
xguidefontsize=40,yguidefontsize=35, color = :teal, label = L"\  (1, 1.5)")
scatter!(time_steps[2:4:end], mQ2[2:4:end], markersize=10, markerstrokewidth=1, color = :lightcoral, label = L"\ (0.1, 0.75)")
scatter!(time_steps[1:4:end], mQ3[1:4:end], markersize=10, markerstrokewidth=1, color = :lightblue,label = L"\ (1, 0.75)")
plot!(xticks=false)
plot!(legend=false)
plot!(widen=false)
vline!([xlims(p1)[2]], lc=:black, lw=2, label = false)
hline!([ylims(p1)[2]], lc=:black, lw=2, label = false)
plot!(grid = false)
ticks_length!(tl=.04)
#= 
curves!(twinx(), x, 2 * cos.(x), yaxis = "Y label 2"; kw...)
 =#

p2 = scatter(time_steps[1:2:end], vQ1[1:2:end], markersize=10, markerstrokewidth=1, color = :teal, seriesalpha = 1, xaxis = L"tω_C", linewidth = 7, label = L"\  (1, 1.5)")
scatter!(time_steps[2:2:end], vQ2[2:2:end], markersize=10, markerstrokewidth=1, color = :lightcoral,seriesalpha = 1, linewidth = 7, label = L"\ (0.1, 0.75)")
scatter!(time_steps[1:2:end], vQ3[1:2:end], markersize=10, markerstrokewidth=1, color = :lightblue,xtickfont = font(24), ytickfont = font(24), xticks = (xticks, xtick_labels),  yticks = (yticks2, ytick_labels2),
seriesalpha = 1, linewidth = 7, label = L"\ (1, 0.75)",xguidefontsize=40,yguidefontsize=35, yaxis = L"\langle\langle Q^2 \rangle\rangle")
#plot!(yaxis=:right)
plot!(widen=false)
vline!([xlims(p2)[2]], lc=:black, lw=2, label = false)
hline!([ylims(p2)[2]], lc=:black, lw=2, label = false)
plot!(grid = false)
plot!(legend = false)
#plot!(legend = :outertopright)
#= plot!(legend = (.6,0.7), grid = false, legendfontsize = 13) #, left_margin = [2mm 0mm], right_margin = [0mm 3mm], bottom_margin = [3mm 3mm])   #,framestyle = :box), legendposition = (1.6, 0.7)
plot!(legendtitle=L"(T,\ \alpha)", legendtitlefontsize=15)
 =#
ticks_length!(tl=.04)

# Create an empty plot for legend only (no axes, no frame)
p3 = plot(framestyle=:none, grid=false, legendfontsize=25)
scatter!([], [], markersize=1, markerstrokewidth=.1, color=:teal, label=L"\ \ \  1.0 \ \  1.50")
scatter!([], [], markersize=1, markerstrokewidth=.1, color=:lightcoral, label=L"\ \ \ 0.1 \ \ 0.75")
scatter!([], [], markersize=1, markerstrokewidth=.1, color=:lightblue, label=L"\ \ \ 1.0 \ \ 0.75")
plot!(legendtitle=L"\ \ \ T\ \ \ \ \alpha", legendtitlefontsize=32)
plot!(foreground_color_legend = nothing, background_color_legend = :white)

#custom_layout = @layout [[a;b] c{0.25w}]  # P1 takes 60% of the column height
#p=plot(p1, p2, p3, layout=custom_layout, size = (900, 800), margin = 5mm)

custom_layout = @layout [[a;b] c{0.25w}]  # P1 takes 60% of the column height
p=plot(p1, p2, p3, layout=custom_layout, size = (900, 800), margin = 5mm)

display("image/png",p)

 