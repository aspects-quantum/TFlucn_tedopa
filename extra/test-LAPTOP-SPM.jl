using DrWatson, Plots, LaTeXStrings
gr()  # Use GR backend, but you can switch to PyPlot if needed

# Define data
a_plus = [0.09495992168234713, 0.33702123217274094, 0.6372444896736239, 0.9231367705835735, 
          1.1615464001371252, 1.347733484238672, 1.4890100635883856, 1.595849570599494,
          1.677332024145704, 1.7402102037323495, 1.790046901376592, 1.8304485200765086, 
          1.8634051422401667, 1.890939939037258, 1.9149037738675725, 1.93495135016066,
          1.9527759838347123, 1.9679265433206998, 1.9808005761639165, 1.9926380385688867, 
          2.002799904210332, 2.012851674669058, 2.0217916794514443, 2.030263953987257, 
          2.038063567343079, 2.044545032604059, 2.0505445534867364, 2.0575423103798114, 
          2.0640886354658194, 2.0696213998269517, 2.0751770882007134]

a_up = [0.09552011750390858, 0.34278645428748716, 0.6586830462975073, 0.9725049773940453,
        1.2483702020746037, 1.477447258343144, 1.6632476355539043, 1.813411162889399,
        1.9355827522857658, 2.035944679745456, 2.1189882271626956, 2.188229279321799,
        2.2468506799215917, 2.2974034397288943, 2.341680621541797, 2.3795126059699805,
        2.4135375490451954, 2.4439780949331515, 2.4708301922974116, 2.4948263635099726,
        2.5164834344953952, 2.5357876912511426, 2.553564175353163, 2.5699792516313194,
        2.5851434156168316, 2.5997227095899067, 2.613347278374249, 2.6258577431914243,
        2.637438606399008, 2.648724372118132, 2.658218253384091]

# Time step duration
tau = 0.04
time_steps = collect(0:length(a_plus)-1) * tau
yticks = collect(range(0,findmax(a_plus)[1]+.2, 5))
xtick_labels = [L"$(x)" for x in time_steps]
ytick_labels = [L"$y" for y in yticks]


# Filter out any potential invalid values
time_steps = filter(x -> isfinite(x) && x >= 0, time_steps)
a_plus = filter(isfinite, a_plus)
a_up = filter(isfinite, a_up)

# Plot with adjusted settings
p1 = plot(time_steps, real(a_plus), label = L"\alpha = 0.5, \vert + \rangle", 
          seriesalpha = 0.7, linewidth = 2, legendfontsize = 16, size=(900, 700))

plot!(time_steps, real(a_up), label = L"\alpha = 0.5, \vert 1 \rangle", 
      seriesalpha = 0.7, linewidth = 2)

# Customize axes labels and ticks
xlabel!(L"t \, (\mathrm{Time\, [s]})", fontsize = 20)
ylabel!(L"\langle Q \rangle \, (\mathrm{Observable})", fontsize = 20)
#= 
# Set tick and label font sizes
plot!(xtickfont = font(16), ytickfont = font(16), guidefont = font(20))
plot!(xticks = (time_steps, xtick_labels), yticks = (yticks, ytick_labels),
      xtickfont = font(16), ytickfont = font(16), guidefont = font(20))
 =#

# Customize plot appearance
plot!(title = "Expectation Value of Q vs Time", titlefontsize = 20,
      legend = :topright, grid = false, framestyle = :box)

# Display and save the plot
display(p1)
file_name_png = "test_1.png"
safesave(file_name_png, p1)
