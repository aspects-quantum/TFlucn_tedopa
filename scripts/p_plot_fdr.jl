
using DrWatson, Plots, LaTeXStrings, PGFPlotsX
#= 
default(
	fontfamily = "Times",
	legendfontsize = 10,
	guidefontsize = 30,
	tickfontsize = 30,
	titlefontsize = 30,
	grid = false,
	linewidth = 10,
	framestyle = :box,
	size = (600, 300),
	dpi = 600,
	minorgrid = false,
	legend = :bottomright,
)
 =#
T_list = [3,5,8,10]

FDR_75 = [4.058477050358806, 2.855138012573611, 2.455431483573631, 2.3563665446356485]
vQ_75 = [44.68176147454178,  44.48321487192148, 46.456515447226685,  46.3026104819018]
mQ_75 = [3.6698299156815906, 3.116011532614108, 2.364987363626919,  1.9650003344051599]

FDR_1 = [9.309102149435843, 5.192141411859455, 3.685828405488081, 3.2847183116593666]
vQ_1 = [122.58251187300293, 91.8283488136699,  72.12150256710174,  62.00233775700252]
mQ_1 = [ 4.3893424559185075, 3.537205231118061,  2.445905459805017, 1.8875998449218716]


@pgf p=Axis({xlabel = L"T", ylabel = L"\frac{\langle\langle Q^2 \rangle\rangle_{\infty}}{T \langle Q \rangle_{\infty}}"},
    PlotInc({ultra_thick, pink, no_markers, opacity = 1}, Coordinates(T_list,FDR_75)),
    LegendEntry(L"\alpha = 0.75"),
    PlotInc({ultra_thick, brown, no_markers}, Coordinates(T_list,FDR_1)),
    LegendEntry(L"\alpha = 1"),
    PlotInc({dashed, no_markers}, Coordinates(T_list,vec(2*ones(length(T_list),1))))
)
pgfsave("fdr.pdf",p; dpi = 600)
#= 
x = T_list

@pgf GroupPlot(
    {
        group_style =
        {
            group_size="2 by 1",
            #= xticklabels_at="edge bottom",
            yticklabels_at="edge left" =#
        },
        #no_markers
    },
    {},
    PlotInc(Coordinates(T_list,FDR_75)),
    LegendEntry(L"\frac{\langle\langle Q^2 \rangle\rangle_{\infty}}{T \langle Q \rangle_{\infty}}, \alpha = 0.75"),
    PlotInc(Coordinates(T_list,FDR_1)),
    LegendEntry(L"\frac{\langle\langle Q^2 \rangle\rangle_{\infty}}{T \langle Q \rangle_{\infty}}, \alpha = 1"),
    {},
    PlotInc(Coordinates(T_list,vQ_75)),
    LegendEntry(L"\langle\langle Q^2 \rangle\rangle_{\infty}, \alpha = 0.75"),
    PlotInc(Coordinates(T_list,vQ_1)),
    LegendEntry(L"\langle\langle Q^2 \rangle\rangle_{\infty}, \alpha = 1"),
    PlotInc(Coordinates(T_list,mQ_75)),
    LegendEntry(L"\langle Q \rangle_{\infty}, \alpha = 0.75"),
    PlotInc(Coordinates(T_list,mQ_1)),
    LegendEntry(L"\langle Q \rangle_{\infty}, \alpha = 1")
) =#



