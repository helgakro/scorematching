####Plotting functions


set_plot_theme <- function(){
  theme_set(theme_bw())
  theme_update(text = element_text(family = "serif", size=12),plot.title = element_text(hjust = 0.5))
}



#' @import cowplot ggplot2
#'
#' @export
plot_grid_2 <- function(p1,p2){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B"),
    hjust = -1,
    nrow = 1
  )

  # extract a legend that is laid out horizontally
  legend_b <- get_legend(
    p2 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )

  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)))
}

#' @import cowplot ggplot2
#'
#' @export
plot_grid_3 <- function(p1,p2,p3){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    p3 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B", "C"),
    hjust = -1,
    nrow = 1
  )

  # extract a legend that is laid out horizontally
  legend_b <- get_legend(
    p2 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )

  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)))
}

plot_grid_3_vertical <- function(p1,p2,p3){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    p3 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B", "C"),
    hjust = -1,
    ncol = 1
  )

  # extract a legend that is laid out horizontally
  legend_b <- get_legend(
    p2 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )

  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)))
}



plot_grid_3_nolegend <- function(p1,p2,p3){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    p3 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B", "C"),
    hjust = -1,
    nrow = 1
  )


  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, ncol = 1, rel_heights = c(1, .1)))
}
# plot(score_par[1,],log_par[1,])
# abline(a=0,b=1)
# plot(score_par[2,],log_par[2,])
# abline(a=0,b=1)


#' @import cowplot ggplot2
#'
#' @export
plot_grid_4 <- function(p1,p2,p3,p4){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    p3 + theme(legend.position="none"),
    p4 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B", "C", "D"),
    hjust = -1,
    nrow = 2
  )

  # extract a legend that is laid out horizontally
  legend_b <- get_legend(
    p4 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )

  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .05)))
}

plot_grid_2_nolegend <- function(p1,p2){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B"),
    hjust = -1,
    nrow = 1
  )


  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, ncol = 1, rel_heights = c(1, .1)))
}
