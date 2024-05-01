## Plot the data on a map, with color indicating values of anomalies

library(fields)
library(maps)
load("anom1962.RData")

color.plot <- function(z, loc, k = 100, zlim = range(z), col = tim.colors(k)[k:1]){  
  breaks <- seq(zlim[1], zlim[2], length = k + 1)
  plot(loc, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  for(i in 1:length(col)){
    index <- z >= breaks[i] & z <= breaks[i+1]
    points(loc[index, 1], loc[index, 2], col = col[i], pch = 16, cex = 0.5)
  }
}

x11(height = 5, width = 8)
par(mar = c(0, 0, 0, 5))
color.plot(z, loc)
map("usa", add = TRUE)
image.plot(legend.only = TRUE, zlim = range(z), col = tim.colors(100)[100:1], legend.mar = 3)
dev.print(device = pdf, file = "anom1962color.pdf", height = 5, width = 8)
dev.print(device = jpeg, file = "anom1962color.jpeg", height = 500, width = 800)

col <- gray(seq(1, 0, length = 100))
color.plot(z, loc, col = col)
map("usa", add = TRUE)
image.plot(legend.only = TRUE, zlim = range(z), col = col, legend.mar = 3)
dev.print(device = pdf, file = "anom1962bw.pdf", height = 5, width = 8)
dev.print(device = jpeg, file = "anom1962bw.jpeg", height = 500, width = 800)

dev.off()



