
recurrence_plot = function(recurrence_matrix){
  # recurrence matrix is the ngCMatrix output from the crqa::crqa function
  
  # # convert the numeric logical matrix into something that R can plot

  # plot_matrix = as(recurrence_matrix, 'matrix')*1 # the *1 converts from logical to numeric
  plot_matrix = t(recurrence_matrix)
  diag(plot_matrix) = 1
  Size = dim(plot_matrix)
  rows = Size[1]
  cols = Size[2]
  oldpar = par(no.readonly = TRUE)
  par(pty = 's')
  par(mgp = c(0,0.5,0))
  par(mar = c(3,1,2,1))
  graphics::image(1:rows,1:cols,plot_matrix,
                  asp = 1,
                  bty = 'n',
                  axes = F, 
                  xlab = '',
                  ylab = '',
                  zlim = c(0,1),
                  col = c(0,1),
                  main = 'Recurrence Plot',
                  useRaster = TRUE)
  
  axis(side = 1, tck = -0.02, xlim = c(0,cols))
  axis(side = 2, tck = -0.02, ylim =c(0,rows))
  axis(side = 3, tck = -0.02, ylim = c(0,cols), labels = FALSE)
  axis(side = 4, tck = -0.02, ylim = c(0,rows), labels = FALSE)
  box(lwd = 1)
  title(ylab = 'Time', line = 1.75)
  title(xlab = 'Time', line = 1.75)
  
  par(oldpar)
}


# library(crqa)
# set.seed(123456)
# n = 100
# test_signal_1 = randi(10, n, 1)
# test_signal_2 = randi(10, n, 1)
# out = crqa(test_signal_1,test_signal_2,
#            delay = 1,
#            embed = 1,
#            rescale = 1,
#            radius = 0.01,
#            normalize = 1,
#            mindiagline = 2,
#            minvertline = 2,
#            tw = 0,
#            whiteline = F,
#            recpt = F)
# 
# recurrence_plot(out$RP)
