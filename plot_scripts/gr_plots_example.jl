using Plots
gr()

kw = (; lab = "", title_loc = :left)
x = π:0.1:2π
#curves([1, 2, 3, 4], [1, 1, 2, 4], title = "Bézier curve")
curves(x, sin.(x), xaxis = "common X label", yaxis = "Y label 1", color = :red, title = "twinx")
p = curves!(twinx(), x, 2 * cos.(x), yaxis = "Y label 2"; kw...)
plot(p)