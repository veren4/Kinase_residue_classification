plot_umap = function(x,
                     labels,
                     # main,
                     colors = c("#ff7f00", "#e377c2", "#17becf"),
                     pad = 0.1, cex = 0.65, pch = 19, add = FALSE, legend.suffix = "",
                     cex.main = 1, cex.legend = 1) {
  
  layout = x
  if (class(x) == "umap") {
    layout = x$layout
  }
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar = c(0.2, 0.7, 1.2, 0.7), ps = 10)   # bottom, left, top, right     0.2, 0.7, 1.2, 0.7
    plot(xylim, xylim, type = "n", axes = F, frame = F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd = 0.25)
  }
  points(layout[,1], layout[,2], col = colors[as.integer(labels)],
         cex = cex, pch = pch)
  
  # mtext(side = 3, main, cex = cex.main)   # title
  
  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  
  # if (add) {
  #   legend.pos = "bottomright"
  #   legend.text = paste(as.character(labels.u), legend.suffix)
  # }
  
  legend(legend.pos, legend = legend.text,
         col = colors[as.integer(labels.u)],
         bty = "n", pch = pch, cex = cex.legend,
         inset = 0.05)
  
}
