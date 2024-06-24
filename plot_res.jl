# Copyright (C) 2024 nikoloz katchakhidze
# Created: 2024-03-21
using PyPlot

function plot_res(fs, xx, fex, fapp, lg1, lg2, ylab)

  _ = figure(figsize=fs)  
  _ = PyPlot.plot(xx,fex, linewidth=3, color="blue") 
  _ = PyPlot.plot(xx,fapp, linewidth=3, color="red", linestyle=(0, (5, 5))) 
  ax = gca()
  axis("tight")
  ax.spines["top"].set_visible(false) # Hide the top edge of the axis
  ax.spines["right"].set_visible(false) # Hide the right edge of the axis
  ax.xaxis.set_ticks_position("bottom") # Set the x-ticks to only the bottom
  ax.yaxis.set_ticks_position("left") # Set the y-ticks to only the left
  ax.legend([lg1, lg2])
  ax.grid(which="both", alpha = 0.1)
  ax.minorticks_on()
  ax.tick_params(which = "minor", bottom = false, left = false)
  xlabel("x")
  ylabel(ylab)
  
end
