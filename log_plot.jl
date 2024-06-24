using PyPlot

function log_plot(fs, x, y, h, x_min, x_max, y_min, y_max, ylab)

  xlab = "number of Galerkin coefficients  n"
  _ = figure(figsize=fs)  
  #_ = PyPlot.plot(x, y, linewidth=3, color="blue")
  _ = PyPlot.plot(x, y, ls = "none", marker="o", ms=5; mec="blue")    
  ax = gca()
  axis("tight")
  ax.spines["top"].set_visible(false) # Hide the top edge of the axis
  ax.spines["right"].set_visible(false) # Hide the right edge of the axis
  ax.xaxis.set_ticks_position("bottom") # Set the x-ticks to only the bottom
  ax.yaxis.set_ticks_position("left") # Set the y-ticks to only the left  
  ax.grid(which="both", alpha = 0.25)
  ax.minorticks_on()
  ax.tick_params(which = "minor", bottom = false, left = false)
  ax.set_xlim(xmin=x_min, xmax=x_max)
  ax.set_ylim(ymin=y_min, ymax=y_max) 
  ax.set_title("h = " * string(h))  
  xscale("log")
  yscale("log")
  xlabel(xlab)
  ylabel("algorithm error  " * ylab)
  
end