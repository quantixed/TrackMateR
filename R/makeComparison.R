#' Make Comparison Plots
#'
#' A series of ggplots to compare between conditions.
#' Called from `compareDatasets()` this function generates plots using summary data of datasets, per condition.
#'
#' @param df data frame called megareport
#' @param msddf data frame of msd averages per dataset
#' @param units character vector of space and time units for axes
#' @param msdplot keyword used to specify msdplot axes, either "loglog", "linlin" (default), "loglin" or "linlog"
#' @param titleStr string used as the title for the report
#' @param subStr string used as the subtitle for the report
#' @return patchwork ggplot
#' @export
makeComparison <- function (df, msddf, units = c("um","s"), msdplot = "linlin", titleStr = "Comparison", subStr = NULL) {
  condition <- dee <- neighbours <- speed <- fd <- width <- NULL

  oldw <- getOption("warn")
  options(warn = -1)

  symlim <- findLog2YAxisLimits(df$alpha)

  # plot alpha comparison
  p_alpha <- ggplot(data = df, aes(x = condition, y = alpha, colour = condition)) +
    geom_boxplot(colour = "grey", outlier.shape = NA) +
    geom_sina(alpha = 0.5, stroke = 0) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey") +
    scale_y_continuous(limits = symlim, trans = "log2") +
    guides(x =  guide_axis(angle = 90)) +
    labs(x = "", y = "Mean alpha") +
    theme_classic() +
    theme(legend.position = "none")

  # plot speed comparison
  p_speed <- ggplot(data = df, aes(x = condition, y = speed, colour = condition)) +
    geom_boxplot(colour = "grey", outlier.shape = NA) +
    geom_sina(alpha = 0.5, stroke = 0) +
    ylim(c(0,NA)) +
    guides(x =  guide_axis(angle = 90)) +
    labs(x = "", y = paste0("Mean speed (",units[1],"/",units[2],")")) +
    theme_classic() +
    theme(legend.position = "none")

  # plot diffusion constant
  p_dee <- ggplot(data = df, aes(x = condition, y = dee, colour = condition)) +
    geom_boxplot(colour = "grey", outlier.shape = NA) +
    geom_sina(alpha = 0.5, stroke = 0) +
    ylim(c(0,NA)) +
    guides(x =  guide_axis(angle = 90)) +
    labs(x = "", y = substitute(paste("Diffusion coefficient (",mm^2,"/",nn,")"), list(mm = units[1], nn = units [2]))) +
    theme_classic() +
    theme(legend.position = "none")

  # plot neighbour density
  p_density <- ggplot(data = df, aes(x = condition, y = neighbours, colour = condition)) +
    geom_boxplot(colour = "grey", outlier.shape = NA) +
    geom_sina(alpha = 0.5, stroke = 0) +
    ylim(c(0,NA)) +
    guides(x =  guide_axis(angle = 90)) +
    labs(x = "", y = paste0("Median track density")) +
    theme_classic() +
    theme(legend.position = "none")

  # plot fractal dimension
  p_fd <- ggplot(data = df, aes(x = condition, y = fd, colour = condition)) +
    geom_boxplot(colour = "grey", outlier.shape = NA) +
    geom_sina(alpha = 0.5, stroke = 0) +
    ylim(c(0,NA)) +
    guides(x =  guide_axis(angle = 90)) +
    labs(x = "", y = paste0("Median fractal dimension")) +
    theme_classic() +
    theme(legend.position = "none")

  # plot width of track summary
  p_width <- ggplot(data = df, aes(x = condition, y = width, colour = condition)) +
    geom_boxplot(colour = "grey", outlier.shape = NA) +
    geom_sina(alpha = 0.5, stroke = 0) +
    ylim(c(0,NA)) +
    guides(x =  guide_axis(angle = 90)) +
    labs(x = "", y = paste0("Median track width")) +
    theme_classic() +
    theme(legend.position = "none")

  # plot msd summary curves altogether
  p_msd <- ggplot(data = msddf, aes(x = t, y = mean, fill = condition)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2) +
    geom_line(aes(colour = condition), size = 1)
  if(msdplot == "linlin") {
    p_msd <- p_msd + ylim(0,NA) + xlim(0,NA)
  } else if(msdplot == "loglin") {
    p_msd <- p_msd + scale_y_log10() + xlim(0,NA)
  } else if(msdplot == "linlog") {
    p_msd <- p_msd + ylim(0,NA) + scale_x_log10()
  } else if(msdplot == "loglog") {
    p_msd <- p_msd + scale_y_log10() + scale_x_log10()
  } else {
    # no scaling applied if any other string is used
  }
  p_msd <- p_msd + labs(x = "Time (s)", y = "MSD") +
    theme_classic() +
    theme(legend.position = "none")

  r_report <- (p_alpha + p_speed + p_dee + p_fd) / (p_width + p_density + p_msd)
  r_report <- r_report + plot_annotation(title = titleStr, subtitle = subStr)

  options(warn = oldw)

  return(r_report)
}
