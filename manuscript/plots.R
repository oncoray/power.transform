get_annotation_settings <- function(ggtheme = NULL) {
  # Import formatting settings from the provided ggtheme.

  # Find the text size for the table. This is based on text sizes in the
  # ggtheme.
  fontsize <- ggtheme$text$size
  fontsize_rel <- 1.0

  # Attempt to base the text size on the general axis.text attribute.
  if (!is.null(ggtheme$axis.text$size)) {
    if (inherits(ggtheme$axis.text$size, "rel")) {
      # Find the relative text size of axis text.
      fontsize_rel <- as.numeric(ggtheme$axis.text$size)
    } else {
      # Set absolute text size.
      fontsize <- ggtheme$axis.text$size
      fontsize_rel <- 1.0
    }
  }

  # Attempt to refine the text size using the axis.text.y attribute in
  # particular.
  if (!is.null(ggtheme$axis.text.y$size)) {
    if (inherits(ggtheme$axis.text.y$size, "rel")) {
      # Set relative text size of axis text
      fontsize_rel <- as.numeric(ggtheme$axis.text.y$size)
    } else {
      # Set absolute text size.
      fontsize <- as.numeric(ggtheme$axis.text.y$size)
      fontsize_rel <- 1.0
    }
  }

  # Update the text size using the magical ggplot2 point size (ggplot2:::.pt).
  geom_text_size <- fontsize * fontsize_rel / 2.845276

  # Obtain lineheight
  lineheight <- ggtheme$text$lineheight
  if (!is.null(ggtheme$axis.text$lineheight)) lineheight <- ggtheme$axis.text$lineheight
  if (!is.null(ggtheme$axis.text.y$lineheight)) lineheight <- ggtheme$axis.text.y$lineheight

  # Obtain family
  fontfamily <- ggtheme$text$family
  if (!is.null(ggtheme$axis.text$family)) fontfamily <- ggtheme$axis.text$family
  if (!is.null(ggtheme$axis.text.y$family)) fontfamily <- ggtheme$axis.text.y$family
  if (!is.null(ggtheme$axis.text.x$family)) fontfamily <- ggtheme$axis.text.x$family

  # Obtain face
  fontface <- ggtheme$text$face
  if (!is.null(ggtheme$axis.text$face)) fontface <- ggtheme$axis.text$face
  if (!is.null(ggtheme$axis.text.y$face)) fontface <- ggtheme$axis.text.y$face
  if (!is.null(ggtheme$axis.text.x$face)) fontface <- ggtheme$axis.text.x$face

  # Obtain colour
  colour <- ggtheme$text$colour
  if (!is.null(ggtheme$axis.text$colour)) colour <- ggtheme$axis.text$colour
  if (!is.null(ggtheme$axis.text.y$colour)) colour <- ggtheme$axis.text.y$colour
  if (!is.null(ggtheme$axis.text.x$colour)) colour <- ggtheme$axis.text.x$colour

  return(list(
    "geom_text_size" = geom_text_size,
    "fontsize" = fontsize,
    "fontsize_rel" = fontsize_rel,
    "colour" = colour,
    "family" = fontfamily,
    "face" = fontface,
    "lineheight" = lineheight
  ))
}



.plot_reduced_normality <- function(plot_theme, manuscript_dir) {
  # Creates the plot that shows how transformation parameter lambda depends on
  # the location for conventional power transformations.

  require(patchwork)
  require(ggplot2)

  # Lambda plot,
  .create_lambda_shift_plot <- function(
      data,
      plot_theme,
      limits,
      guide = FALSE,
      strip_y_axis = TRUE) {
    # Prevent warnings due to non-standard evaluation.
    d <- lambda <- method <- version <- NULL

    annotation_settings <- get_annotation_settings(plot_theme)

    p <- ggplot2::ggplot(
      data = data,
      mapping = ggplot2::aes(
        x = d,
        y = lambda,
        colour = method
      )
    )
    p <- p + plot_theme
    p <- p + ggplot2::geom_point()
    p <- p + ggplot2::geom_hline(
      yintercept = 1.0,
      linetype = "longdash",
      colour = "grey40"
    )
    p <- p + ggplot2::annotate(
      geom = "text",
      x = 6,
      y = 1,
      label = "expected",
      colour = annotation_settings$colour,
      family = annotation_settings$family,
      fontface = annotation_settings$face,
      size = annotation_settings$geom_text_size,
      vjust = -1.0,
      hjust = 1.0
    )
    p <- p + ggplot2::scale_x_continuous(
      name = latex2exp::TeX("$\\mu$"),
      labels = scales::math_format()
    )
    p <- p + ggplot2::scale_y_continuous(
      name = latex2exp::TeX("$\\lambda$"),
      limits = limits
    )

    if (guide) {
      p <- p + paletteer::scale_color_paletteer_d(
        palette = "ggthemes::Tableau_10",
        drop = FALSE
      )
      p <- p + guides(colour = ggplot2::guide_legend(title = "method"))
    } else {
      p <- p + paletteer::scale_color_paletteer_d(
        palette = "ggthemes::Tableau_10",
        drop = FALSE,
        guide = "none"
      )
    }

    if (strip_y_axis) {
      p <- p + ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
    }

    return(p)
  }

  # Prevent warnings due to non-standard evaluation.
  estimation_method <- distribution <- method <- version <- NULL

  # Process / read data.
  data <- .get_shifted_distribution_data(manuscript_dir = manuscript_dir)
  data <- data[estimation_method == "MLE" & version == "original" & distribution == "normal"]

  # Normal distribution --------------------------------------------------------

  # Box-Cox
  p_bc_normal <- .create_lambda_shift_plot(
    data = data[method == "Box-Cox"],
    plot_theme = plot_theme,
    limits = c(-5.0, 35.0),
    strip_y_axis = FALSE
  )

  # Yeo-Johnson
  p_yj_normal <- .create_lambda_shift_plot(
    data = data[method == "Yeo-Johnson"],
    plot_theme = plot_theme,
    limits = c(-5.0, 35.0),
    guide = TRUE
  )

  # Patch all the plots together.
  p <- p_bc_normal + p_yj_normal + patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  )

  return(p)
}



.plot_shifted_distributions <- function(plot_theme, manuscript_dir) {
  # Creates the plot for lambda values under location shift for MLE

  require(patchwork)
  require(ggplot2)

  .create_density_plot <- function(
      x,
      plot_theme,
      limits) {
    p <- ggplot2::ggplot(
      data = data.table::data.table("x" = x),
      mapping = ggplot2::aes(x = x)
    )
    p <- p + plot_theme
    p <- p + ggplot2::geom_density()
    p <- p + ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )
    p <- p + ggplot2::xlim(limits)

    return(p)
  }

  # Lambda plot,
  .create_lambda_shift_plot <- function(
      data,
      plot_theme,
      limits,
      guide = FALSE,
      strip_x_axis = TRUE,
      strip_y_label = TRUE) {
    # Prevent warnings due to non-standard evaluation.
    d <- lambda <- method <- version <- NULL

    p <- ggplot2::ggplot(
      data = data,
      mapping = ggplot2::aes(
        x = d,
        y = lambda,
        colour = method,
        shape = version
      )
    )
    p <- p + plot_theme
    p <- p + ggplot2::geom_point()
    p <- p + ggplot2::scale_x_continuous(
      name = latex2exp::TeX("$d$"),
      labels = scales::math_format()
    )
    p <- p + ggplot2::scale_y_continuous(
      name = latex2exp::TeX("$\\lambda$"),
      limits = limits
    )

    if (guide) {
      p <- p + paletteer::scale_color_paletteer_d(
        palette = "ggthemes::Tableau_10",
        drop = FALSE
      )
      p <- p + guides(colour = ggplot2::guide_legend(title = "method"))
    } else {
      p <- p + paletteer::scale_color_paletteer_d(
        palette = "ggthemes::Tableau_10",
        drop = FALSE,
        guide = "none"
      )
      p <- p + ggplot2::scale_shape_discrete(guide = "none")
    }

    if (strip_x_axis) {
      p <- p + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    }

    if (strip_y_label) {
      p <- p + ggplot2::theme(
        axis.title.y = ggplot2::element_blank()
      )
    }

    return(p)
  }

  # Prevent warnings due to non-standard evaluation.
  estimation_method <- distribution <- method <- version <- NULL

  # Set seed.
  set.seed(19L)

  # Normal distribution.
  x_normal <- power.transform::ragn(
    10000L,
    location = 0,
    scale = 1 / sqrt(2),
    alpha = 0.5,
    beta = 2
  )

  # Right skewed data
  x_right_skewed <- power.transform::ragn(
    10000L,
    location = 0,
    scale = 1 / sqrt(2),
    alpha = 0.2,
    beta = 2
  )

  # Left skewed data
  x_left_skewed <- power.transform::ragn(
    10000L,
    location = 0,
    scale = 1 / sqrt(2),
    alpha = 0.8,
    beta = 2
  )

  # Process / read data.
  data <- .get_shifted_distribution_data(manuscript_dir = manuscript_dir)
  data <- data[estimation_method == "MLE"]

  # Normal distribution --------------------------------------------------------

  # Density plot
  p_dens_normal <- .create_density_plot(
    x = x_normal,
    plot_theme = plot_theme,
    limits = c(-6.0, 6.0)
  )
  p_dens_normal <- p_dens_normal + ggplot2::facet_grid(cols = vars("normal distribution"))

  # Box-Cox
  p_bc_normal <- .create_lambda_shift_plot(
    data = data[distribution == "normal" & method == "Box-Cox"],
    plot_theme = plot_theme,
    limits = c(-5.0, 35.0),
    strip_y_label = FALSE
  )

  # Yeo-Johnson
  p_yj_normal <- .create_lambda_shift_plot(
    data = data[distribution == "normal" & method == "Yeo-Johnson"],
    plot_theme = plot_theme,
    limits = c(-5.0, 35.0),
    strip_y_label = FALSE,
    strip_x_axis = FALSE
  )

  # Right skewed distribution --------------------------------------------------

  # Density
  p_dens_right <- .create_density_plot(
    x = x_right_skewed,
    plot_theme = plot_theme,
    limits = c(-3.0, 10.0)
  )
  p_dens_right <- p_dens_right + ggplot2::facet_grid(cols = vars("right-skewed distribution"))

  # Box-Cox-original
  p_bc_right <- .create_lambda_shift_plot(
    data = data[distribution == "right-skewed" & method == "Box-Cox"],
    plot_theme = plot_theme,
    limits = c(-5.0, 1.0)
  )

  # Yeo-Johnson
  p_yj_right <- .create_lambda_shift_plot(
    data = data[distribution == "right-skewed" & method == "Yeo-Johnson"],
    plot_theme = plot_theme,
    limits = c(-5.0, 1.0),
    strip_x_axis = FALSE
  )

  # Left skewed distribution ---------------------------------------------------

  # Density
  p_dens_left <- .create_density_plot(
    x = x_left_skewed,
    plot_theme = plot_theme,
    limits = c(-10.0, 3.0)
  )
  p_dens_left <- p_dens_left + ggplot2::facet_grid(cols = vars("left-skewed distribution"))

  # Box-Cox
  p_bc_left <- .create_lambda_shift_plot(
    data = data[distribution == "left-skewed" & method == "Box-Cox"],
    plot_theme = plot_theme,
    limits = c(0.0, 65.0),
    guide = TRUE
  )

  # Yeo-Johnson
  p_yj_left <- .create_lambda_shift_plot(
    data = data[distribution == "left-skewed" & method == "Yeo-Johnson"],
    plot_theme = plot_theme,
    limits = c(0.0, 65.0),
    strip_x_axis = FALSE
  )

  # Patch all the plots together.
  p <- p_dens_normal + p_dens_right + p_dens_left +
    p_bc_normal + p_bc_right + p_bc_left +
    p_yj_normal + p_yj_right + p_yj_left +
    patchwork::plot_layout(
      ncol = 3,
      heights = c(0.5, 1.0, 1.0),
      guides = "collect"
    )

  return(p)
}



.plot_shifted_distributions_other_criteria <- function(plot_theme, manuscript_dir) {
  # Creates the plot for lambda values under location shift for all optimisation
  # criteria.

  require(patchwork)
  require(ggplot2)

  # Lambda plot,
  .create_lambda_shift_plot <- function(
      data,
      plot_theme,
      limits,
      guide = FALSE,
      strip_x_axis = TRUE,
      strip_y_label = TRUE) {
    # Prevent warnings due to non-standard evaluation.
    d <- lambda <- estimation_method <- NULL

    p <- ggplot2::ggplot(
      data = data,
      mapping = ggplot2::aes(
        x = d,
        y = lambda,
        colour = estimation_method
      )
    )
    p <- p + plot_theme
    p <- p + ggplot2::geom_point()
    p <- p + ggplot2::scale_x_continuous(
      name = latex2exp::TeX("$d$"),
      labels = scales::math_format()
    )
    p <- p + ggplot2::scale_y_continuous(
      name = latex2exp::TeX("$\\lambda$"),
      limits = limits
    )

    if (guide) {
      p <- p + paletteer::scale_color_paletteer_d(
        palette = "ggthemes::Tableau_10",
        drop = FALSE
      )
      p <- p + guides(colour = ggplot2::guide_legend(title = "optimisation criterion"))
    } else {
      p <- p + paletteer::scale_color_paletteer_d(
        palette = "ggthemes::Tableau_10",
        drop = FALSE,
        guide = "none"
      )
    }

    if (strip_x_axis) {
      p <- p + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    }

    if (strip_y_label) {
      p <- p + ggplot2::theme(
        axis.title.y = ggplot2::element_blank()
      )
    }

    return(p)
  }

  # Prevent warnings due to non-standard evaluation.
  distribution <- method <- version <- NULL

  # Process / read data.
  data <- .get_shifted_distribution_data(manuscript_dir = manuscript_dir)

  # Normal distribution --------------------------------------------------------

  # Box-Cox-original
  p_bc_normal_orig <- .create_lambda_shift_plot(
    data = data[distribution == "normal" & method == "Box-Cox" & version == "original"],
    plot_theme = plot_theme,
    limits = c(-5.0, 35.0),
    strip_y_label = FALSE
  )
  p_bc_normal_orig <- p_bc_normal_orig + ggplot2::facet_grid(cols = vars("normal distribution"))


  # Box-Cox-shift-sensitive
  p_bc_normal_ss <- .create_lambda_shift_plot(
    data = data[distribution == "normal" & method == "Box-Cox" & version == "shift-sensitive"],
    plot_theme = plot_theme,
    limits = c(0.5, 1.5),
    strip_y_label = FALSE
  )

  # Yeo-Johnson-original
  p_yj_normal_orig <- .create_lambda_shift_plot(
    data = data[distribution == "normal" & method == "Yeo-Johnson" & version == "original"],
    plot_theme = plot_theme,
    limits = c(-5.0, 35.0),
    strip_y_label = FALSE
  )

  # Yeo-Johnson-shift-sensitive
  p_yj_normal_ss <- .create_lambda_shift_plot(
    data = data[distribution == "normal" & method == "Yeo-Johnson" & version == "shift-sensitive"],
    plot_theme = plot_theme,
    limits = c(0.5, 1.5),
    strip_x_axis = FALSE,
    strip_y_label = FALSE
  )

  # Right skewed distribution --------------------------------------------------

  # Box-Cox-original
  p_bc_right_orig <- .create_lambda_shift_plot(
    data = data[distribution == "right-skewed" & method == "Box-Cox" & version == "original"],
    plot_theme = plot_theme,
    limits = c(-10.0, 1.0)
  )
  p_bc_right_orig <- p_bc_right_orig + ggplot2::facet_grid(cols = vars("right-skewed distribution"))

  # Box-Cox-shift-sensitive
  p_bc_right_ss <- .create_lambda_shift_plot(
    data = data[distribution == "right-skewed" & method == "Box-Cox" & version == "shift-sensitive"],
    plot_theme = plot_theme,
    limits = c(-0.5, 0.5)
  )

  # Yeo-Johnson-original
  p_yj_right_orig <- .create_lambda_shift_plot(
    data = data[distribution == "right-skewed" & method == "Yeo-Johnson" & version == "original"],
    plot_theme = plot_theme,
    limits = c(-10.0, 1.0)
  )

  # Yeo-Johnson-shift-sensitive
  p_yj_right_ss <- .create_lambda_shift_plot(
    data = data[distribution == "right-skewed" & method == "Yeo-Johnson" & version == "shift-sensitive"],
    plot_theme = plot_theme,
    limits = c(0.0, 1.0),
    strip_x_axis = FALSE
  )

  # Left skewed distribution ---------------------------------------------------

  # Box-Cox-original
  p_bc_left_orig <- .create_lambda_shift_plot(
    data = data[distribution == "left-skewed" & method == "Box-Cox" & version == "original"],
    plot_theme = plot_theme,
    limits = c(0.0, 65.0),
    guide = TRUE
  )
  p_bc_left_orig <- p_bc_left_orig + ggplot2::facet_grid(cols = vars("left-skewed distribution"), rows = vars("Box-Cox"))

  # Box-Cox-shift-sensitive
  p_bc_left_ss <- .create_lambda_shift_plot(
    data = data[distribution == "left-skewed" & method == "Box-Cox" & version == "shift-sensitive"],
    plot_theme = plot_theme,
    limits = c(3.5, 4.5)
  )
  p_bc_left_ss <- p_bc_left_ss + ggplot2::facet_grid(rows = vars("shift-sens. Box-Cox"))

  # Yeo-Johnson-original
  p_yj_left_orig <- .create_lambda_shift_plot(
    data = data[distribution == "left-skewed" & method == "Yeo-Johnson" & version == "original"],
    plot_theme = plot_theme,
    limits = c(0.0, 65.0)
  )
  p_yj_left_orig <- p_yj_left_orig + ggplot2::facet_grid(rows = vars("Yeo-Johnson"))

  # Yeo-Johnson-shift-sensitive
  p_yj_left_ss <- .create_lambda_shift_plot(
    data = data[distribution == "left-skewed" & method == "Yeo-Johnson" & version == "shift-sensitive"],
    plot_theme = plot_theme,
    limits = c(1.0, 2.0),
    strip_x_axis = FALSE
  )
  p_yj_left_ss <- p_yj_left_ss + ggplot2::facet_grid(rows = vars("shift-sens. Yeo-Johnson"))

  # Patch all the plots together.
  p <- p_bc_normal_orig + p_bc_right_orig + p_bc_left_orig +
    p_bc_normal_ss + p_bc_right_ss + p_bc_left_ss +
    p_yj_normal_orig + p_yj_right_orig + p_yj_left_orig +
    p_yj_normal_ss + p_yj_right_ss + p_yj_left_ss +
    patchwork::plot_layout(ncol = 3, guides = "collect")

  return(p)
}



.plot_weighting_functions <- function(plot_theme) {
  # Non-standard evaluation
  weight <- NULL

  # Step function --------------------------------------------------------------
  data_step <- data.table::data.table("x" = seq(from = -1.0, to = 1.0, by = 0.001))
  data_step[, "weight" := power.transform:::..weight_function_external(x, "step", k1 = 0.60)]

  p_step <- ggplot2::ggplot(
    data = data_step,
    mapping = ggplot2::aes(
      x = x,
      y = weight
    )
  )
  p_step <- p_step + plot_theme
  p_step <- p_step + ggplot2::geom_line()
  p_step <- p_step + ggplot2::facet_grid(cols = vars("step"))

  # Triangle function ----------------------------------------------------------
  data_triangle <- data.table::data.table("x" = seq(from = -1.0, to = 1.0, by = 0.001))
  data_triangle[, "weight" := power.transform:::..weight_function_external(x, "triangle", k1 = 0.30, k2 = 0.90)]

  p_triangle <- ggplot2::ggplot(
    data = data_triangle,
    mapping = ggplot2::aes(
      x = x,
      y = weight
    )
  )
  p_triangle <- p_triangle + plot_theme
  p_triangle <- p_triangle + ggplot2::geom_line()
  p_triangle <- p_triangle + ggplot2::facet_grid(cols = vars("triangle"))
  p_triangle <- p_triangle + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )

  # Cosine function ------------------------------------------------------------
  data_cosine <- data.table::data.table(x = seq(from = -1.0, to = 1.0, by = 0.001))
  data_cosine[, "weight" := power.transform:::..weight_function_external(x, "cosine", k1 = 0.30, k2 = 0.90)]

  p_cosine <- ggplot2::ggplot(
    data = data_cosine,
    mapping = ggplot2::aes(
      x = x,
      y = weight
    )
  )
  p_cosine <- p_cosine + plot_theme
  p_cosine <- p_cosine + ggplot2::geom_line()
  p_cosine <- p_cosine + ggplot2::facet_grid(cols = vars("tapered cosine"))
  p_cosine <- p_cosine + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )

  p <- p_step + p_triangle + p_cosine + patchwork::plot_layout(ncol = 3)

  return(p)
}



.plot_optimised_weighting_function_parameters <- function(plot_theme, manuscript_dir) {
  # Prevent warnings due to non-standard evaluation.
  estimation_method <- weight_method <- method <- k <- n <- NULL
  lambda <- target_lambda <- lambda_error <- i.lambda_error <- non_robust_lambda_error <- NULL
  better_than_non_robust <- pass_rate <- median_error_round <- NULL

  # This find the optimal weighting parameter settings that are used within the
  # package. The output is not used here, but its presence is required for
  # processing.
  two_sided_function_parameters <- .get_optimised_weighting_function_parameters(
    manuscript_dir = manuscript_dir,
    side = "both"
  )

  # This finds lambda values at the optimised weighting parameter settings.
  data <- .get_transformation_parameters(manuscript_dir = manuscript_dir)

  # Select only MLE here.
  data <- data[estimation_method == "MLE"]

  # Compute difference to target.
  data <- data[, "lambda_error" := abs(lambda - target_lambda)]

  # Join with error for non-robust values.
  non_robust_data <- data[weight_method == "non-robust", mget(c("method", "estimation_method", "k", "ii", "lambda_error"))]
  data[
    non_robust_data,
    on = c("method", "estimation_method", "k", "ii"),
    "non_robust_lambda_error" := i.lambda_error
  ]

  data[, "better_than_non_robust" := lambda_error <= non_robust_lambda_error]
  data[k == 0.0, "better_than_non_robust" := NA_integer_]
  data[weight_method == "non-robust", "better_than_non_robust" := NA_integer_]

  pass_rate_data <- data[, list(
    "pass_rate" = sum(better_than_non_robust, na.rm = TRUE),
    "n" = sum(!is.na(better_than_non_robust)),
    "median_error" = stats::median(lambda_error)
  ),
  by = c("method", "estimation_method", "weight_method")
  ]
  pass_rate_data[, "pass_rate" := round(pass_rate / n * 100.0, 1)]
  pass_rate_data[, "median_error_round" := format(pass_rate_data$median_error, digits = 1)]

  annotation_settings <- get_annotation_settings(plot_theme)

  p <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
      x = weight_method,
      y = lambda_error,
      fill = method
    )
  )
  p <- p + plot_theme
  p <- p + ggplot2::geom_violin(draw_quantiles = 0.5)
  p <- p + ggplot2::facet_wrap("method", nrow = 2)
  p <- p + paletteer::scale_fill_paletteer_d(
    palette = "ggthemes::Tableau_10",
    drop = FALSE
  )
  p <- p + geom_text(
    data = pass_rate_data,
    mapping = ggplot2::aes(
      x = weight_method,
      label = median_error_round,
      y = 10.0
    ),
    colour = annotation_settings$colour,
    family = annotation_settings$family,
    fontface = annotation_settings$face,
    size = annotation_settings$geom_text_size,
    vjust = 0.5
  )
  p <- p + ggplot2::scale_y_sqrt(name = latex2exp::TeX("$| \\hat{\\lambda}^r - \\hat{\\lambda}_{0}|$"))
  p <- p + ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
  p <- p + ggplot2::xlab("weighting method")

  return(p)
}



.plot_residuals <- function(plot_theme, manuscript_dir) {
  # Prevent warnings due to non-standard evaluation.
  has_outliers <- residual <- threshold <- residual_error <- method <- p <- NULL

  if (!file.exists(file.path(manuscript_dir, "residual_error_plot_main_manuscript.RDS"))) {
    # Get data
    data <- .get_goodness_of_fit_data(manuscript_dir = manuscript_dir)

    outlier_free_data <- data[
      has_outliers == FALSE,
      list("residual_error" = stats::quantile(residual, probs = 0.99)),
      by = c("p", "method")
    ]

    data <- data[, list(
      "residual_50" = stats::quantile(residual, probs = 0.50),
      "residual_90" = stats::quantile(residual, probs = 0.90),
      "residual_95" = stats::quantile(residual, probs = 0.95),
      "residual_99" = stats::quantile(residual, probs = 0.99),
      "residual_max" = max(residual)
    ),
    by = c("p", "method")
    ]

    data <- data.table::melt(
      data = data,
      id.vars = c("p", "method"),
      variable.name = "threshold",
      value.name = "residual_error"
    )

    data$threshold <- factor(
      x = data$threshold,
      levels = c("residual_50", "residual_90", "residual_95", "residual_99", "residual_max"),
      labels = c("50 %", "90 %", "95 %", "99 %", "100 %")
    )

    saveRDS(
      data,
      file.path(manuscript_dir, "residual_error_plot_main_manuscript.RDS")
    )
  } else {
    data <- readRDS(file.path(manuscript_dir, "residual_error_plot_main_manuscript.RDS"))

    outlier_free_data <- data[
      has_outliers == FALSE,
      list("residual_error" = stats::quantile(residual, probs = 0.99)),
      by = c("p", "method")
    ]
  }

  p_bc <- ggplot2::ggplot(
    data = data[method == "Box-Cox"],
    mapping = ggplot2::aes(
      x = p,
      y = residual_error,
      colour = threshold
    )
  )
  p_bc <- p_bc + plot_theme
  p_bc <- p_bc + ggplot2::geom_line()
  p_bc <- p_bc + ggplot2::geom_line(
    data = outlier_free_data[method == "Box-Cox"],
    mapping = ggplot2::aes(
      x = p,
      y = residual_error
    ),
    colour = "gray40",
    linetype = "dashed"
  )
  p_bc <- p_bc + ggplot2::scale_colour_discrete(
    name = "percentile (Box-Cox)",
    type = c("50 %" = "#ABC6E2", "90 %" = "#779EC6", "95 %" = "#4E79A7", "99 %" = "#346394", "100 %" = "#1D4D7E")
  )
  p_bc <- p_bc + ggplot2::xlab("empirical probability")
  p_bc <- p_bc + ggplot2::ylab("absolute residual error")
  p_bc <- p_bc + ggplot2::coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )

  p_yj <- ggplot2::ggplot(
    data = data[method == "Yeo-Johnson"],
    mapping = ggplot2::aes(
      x = p,
      y = residual_error,
      colour = threshold
    )
  )
  p_yj <- p_yj + plot_theme
  p_yj <- p_yj + ggplot2::geom_line()
  p_yj <- p_yj + ggplot2::geom_line(
    data = outlier_free_data[method == "Yeo-Johnson"],
    mapping = ggplot2::aes(
      x = p,
      y = residual_error
    ),
    colour = "gray40",
    linetype = "dashed"
  )
  p_yj <- p_yj + ggplot2::scale_colour_discrete(
    name = "percentile (Yeo-Johnson)",
    type = c("50 %" = "#FFBD7D", "90 %" = "#FFA954", "95 %" = "#F28E2B", "99 %" = "#CD6B0B", "100 %" = "#A25000")
  )
  p_yj <- p_yj + ggplot2::xlab("empirical probability")
  p_yj <- p_yj + ggplot2::ylab("absolute residual error")
  p_yj <- p_yj + ggplot2::coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )
  p_yj <- p_yj + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )

  p <- p_bc + p_yj + patchwork::plot_layout(ncol = 2, guides = "collect")

  return(p)
}



.plot_type_1_error_rate <- function(plot_theme, manuscript_dir) {
  # Prevent warnings due to non-standard evaluation.
  residual <- rejected <- p <- n <- mare <- method <- NULL

  if (!file.exists(file.path(manuscript_dir, "type_1_error_rate_plot_main_manuscript.RDS"))) {
    # Get data
    data <- .get_goodness_of_fit_data(manuscript_dir = manuscript_dir)

    central_width <- c(0.60, 0.70, 0.80, 0.90, 0.95, 1.00)

    mare_data <- list()

    for (ii in seq_along(central_width)) {
      p_lower <- 0.50 - central_width[ii] / 2
      p_upper <- 0.50 + central_width[ii] / 2

      x <- data.table::copy(data)

      # Clip empirical probabilities to centre.
      x <- x[p >= p_lower & p <= p_upper]

      # Compute mean absolute residual error per feature.
      x <- x[, list("mare" = mean(residual)), by = c("distribution_id", "outlier_id", "method")]
      x <- x[, list("n" = .N), by = c("mare", "method")][order(mare, method)]
      x[, "rejected" := 1.0 - cumsum(n) / sum(n), by = "method"]
      x <- rbind(
        data.table::data.table(
          "mare" = c(0.0, 0.0),
          "method" = c("Box-Cox", "Yeo-Johnson"),
          "n" = 0L,
          "rejected" = 1.0
        ),
        x
      )

      x[, "kappa" := central_width[ii]]

      mare_data[[ii]] <- x
    }

    data <- data.table::rbindlist(mare_data)

    data$kappa <- factor(
      x = data$kappa,
      levels = central_width,
      labels = c("60 %", "70 %", "80 %", "90 %", "95 %", "100 %")
    )

    saveRDS(
      data,
      file.path(manuscript_dir, "type_1_error_rate_plot_main_manuscript.RDS")
    )
  } else {
    data <- readRDS(file.path(manuscript_dir, "type_1_error_rate_plot_main_manuscript.RDS"))
  }

  p_bc <- ggplot2::ggplot(
    data = data[method == "Box-Cox"],
    mapping = ggplot2::aes(
      x = mare,
      y = rejected,
      colour = kappa
    )
  )
  p_bc <- p_bc + plot_theme
  p_bc <- p_bc + ggplot2::geom_step()
  p_bc <- p_bc + ggplot2::scale_colour_discrete(
    name = "central portion κ\n(Box-Cox)",
    type = c(
      "60 %" = "#bacbde",
      "70 %" = "#98b2cd",
      "80 %" = "#7598bd",
      "90 %" = "#537dac",
      "95 %" = "#42648a",
      "100 %" = "#324b67"
    )
  )
  p_bc <- p_bc + ggplot2::xlab("test statistic")
  p_bc <- p_bc + ggplot2::ylab("type I error rate")
  p_bc <- p_bc + ggplot2::coord_cartesian(xlim = c(0, 0.25))

  p_yj <- ggplot2::ggplot(
    data = data[method == "Yeo-Johnson"],
    mapping = ggplot2::aes(
      x = mare,
      y = rejected,
      colour = kappa
    )
  )
  p_yj <- p_yj + plot_theme
  p_yj <- p_yj + ggplot2::geom_step()
  p_yj <- p_yj + ggplot2::scale_colour_discrete(
    name = "central portion κ\n(Yeo-Johnson)",
    type = c(
      "60 %" = "#f9c59f",
      "70 %" = "#f6a76f",
      "80 %" = "#f38a3f",
      "90 %" = "#f06d0f",
      "95 %" = "#c0570c",
      "100 %" = "#904109"
    )
  )
  p_yj <- p_yj + ggplot2::xlab("test statistic")
  p_yj <- p_yj + ggplot2::coord_cartesian(xlim = c(0, 0.25))
  p_yj <- p_yj + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )

  p <- p_bc + p_yj + patchwork::plot_layout(ncol = 2, guides = "collect")

  return(p)
}



.plot_type_1_error_rate_appendix <- function(plot_theme, manuscript_dir) {
  # Prevent warnings due to non-standard evaluation.
  residual <- rejected <- p <- n <- mare <- method <- NULL

  # Get data
  data <- .get_test_statistics_data_appendix(manuscript_dir = manuscript_dir)

  central_width <- c(0.60, 0.70, 0.80, 0.90, 0.95, 1.00)

  data$kappa <- factor(
    x = data$kappa,
    levels = central_width,
    labels = c("60 %", "70 %", "80 %", "90 %", "95 %", "100 %")
  )

  data <- data[kappa == "80 %"]

  plot_start_list <- list()
  ii <- 1L
  for (method in levels(data$method)) {
    for (method_tag in levels(data$method_tag)) {
      for (kappa in levels(data$kappa)) {
        plot_start_list[[ii]] <- data.table::data.table(
          mare = 0.0,
          method = method,
          method_tag = method_tag,
          n = 0L,
          rejected = 1.0,
          kappa = kappa
        )

        ii <- ii + 1L
      }
    }
  }

  data <- data.table::rbindlist(c(list(data), plot_start_list))

  p_bc <- ggplot2::ggplot(
    data = data[method == "Box-Cox"],
    mapping = ggplot2::aes(
      x = mare,
      y = rejected,
      colour = method_tag
    )
  )
  p_bc <- p_bc + plot_theme
  p_bc <- p_bc + ggplot2::geom_step()
  p_bc <- p_bc + ggplot2::scale_colour_discrete(
    name = "power transform variant\n(Box-Cox)",
    type = c(
      "robust, shift sens. (MLE)" = "#bacbde",
      "conventional (MLE)" = "#98b2cd",
      "robust, shift sens. (C-vM)" = "#42648a",
      "conventional (C-vM)" = "#324b67"
    )
  )
  p_bc <- p_bc + ggplot2::xlab("test statistic")
  p_bc <- p_bc + ggplot2::ylab("type I error rate")
  p_bc <- p_bc + ggplot2::coord_cartesian(xlim = c(0, 0.15))

  p_yj <- ggplot2::ggplot(
    data = data[method == "Yeo-Johnson"],
    mapping = ggplot2::aes(
      x = mare,
      y = rejected,
      colour = method_tag
    )
  )
  p_yj <- p_yj + plot_theme
  p_yj <- p_yj + ggplot2::geom_step()
  p_yj <- p_yj + ggplot2::scale_colour_discrete(
    name = "power transform variant\n(Yeo-Johnson)",
    type = c(
      "robust, shift sens. (MLE)" = "#f9c59f",
      "conventional (MLE)" = "#f6a76f",
      "robust, shift sens. (C-vM)" = "#c0570c",
      "conventional (C-vM)" = "#904109"
    )
  )
  p_yj <- p_yj + ggplot2::xlab("test statistic")
  p_yj <- p_yj + ggplot2::coord_cartesian(xlim = c(0, 0.15))
  p_yj <- p_yj + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )

  p <- p_bc + p_yj + patchwork::plot_layout(ncol = 2, guides = "collect")

  return(p)
}


.plot_bimodal_distribution_test <- function(plot_theme, manuscript_dir) {

  # Lambda plot,
  .create_density_and_qq_plot <- function(
    location_shift = 0.0,
    plot_theme,
    limits = c(-7.0, 7.0),
    drop_qq_y_axis = TRUE) {
    # Prevent warnings due to non-standard evaluation.
    distribution <- NULL

    set.seed(9L)

    # Normal distribution.
    x_1 <- power.transform::ragn(
      10000L,
      location = -location_shift / 2,
      scale = 1 / sqrt(2),
      alpha = 0.5,
      beta = 2
    )

    # Set data.
    x_2 <- x_1 + location_shift
    data <- data.table::rbindlist(list(
      data.table::data.table(x = x_1, distribution = 1L),
      data.table::data.table(x = x_2, distribution = 2L)
    ))

    data$distribution <- factor(data$distribution, levels = c("1", "2"))

    # Compute p-value without changing the distribution.
    transformer <- suppressWarnings(power.transform::find_transformation_parameters(
      x = data$x,
      method = "none"
    ))

    p_value <- suppressMessages(power.transform::assess_transformation(
      x = data$x,
      transformer = transformer
    ))

    if(p_value < 1E-3){
      p_value <- paste0("p < 0.001")
    } else {
      p_value <- paste0("p = ", sprintf(p_value, fmt = "%#.3f"))
    }

    annotation_settings <- get_annotation_settings(plot_theme)

    p <- ggplot2::ggplot(
      data = data)
    p <- p + plot_theme
    p <- p + ggplot2::geom_density(
      mapping = ggplot2::aes(
        x = x,
        y = ggplot2::after_stat(scaled),
        fill = distribution,
        colour = distribution),
      alpha = 0.2)
    p <- p + ggplot2::geom_density(
      mapping = ggplot2::aes(
        x = x,
        y = ggplot2::after_stat(scaled)),
      colour = "grey10")

    p <- p + ggplot2::geom_segment(
      x = -location_shift / 2.0,
      xend = location_shift / 2.0,
      y = 0.10,
      yend = 0.10,
      colour = "grey40")
    p <- p + ggplot2::annotate(
      geom = "text",
      x = 0.0,
      y = 0.1,
      label = p_value,
      colour = "grey40",
      family = annotation_settings$family,
      fontface = annotation_settings$face,
      size = annotation_settings$geom_text_size,
      vjust = 2.0,
      hjust = 0.5
    )

    p <- p + ggplot2::coord_cartesian(xlim = limits)

    p <- p + ggplot2::facet_grid(cols = vars(paste0("d = ", location_shift)))

    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )

    p <- p + paletteer::scale_color_paletteer_d(
      palette = "ggthemes::Tableau_10",
      drop = FALSE,
      guide = "none"
    )

    p <- p + paletteer::scale_fill_paletteer_d(
      palette = "ggthemes::Tableau_10",
      drop = FALSE,
      guide = "none"
    )

    p_qq <- power.transform::plot_qq_plot(
      x = data$x,
      transformer = transformer,
      show_original = FALSE,
      ggtheme = plot_theme)

    p_qq <- p_qq + ggplot2::xlab("Expected quantile")
    p_qq <- p_qq + ggplot2::ylab("Observed quantile")

    p_qq <- p_qq + ggplot2::coord_cartesian(
      xlim = c(-4.0, 4.0),
      ylim = c(-4.0, 4.0)
    )

    if (drop_qq_y_axis) {
      p_qq <- p_qq + ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      )
    }

    return(list(
      "density" = p,
      "qq" = p_qq))
  }

  p_1 <- .create_density_and_qq_plot(location_shift = 1.0, plot_theme = plot_theme, drop_qq_y_axis = FALSE)
  p_2 <- .create_density_and_qq_plot(location_shift = 2.0, plot_theme = plot_theme)
  p_3 <- .create_density_and_qq_plot(location_shift = 3.0, plot_theme = plot_theme)
  p_4 <- .create_density_and_qq_plot(location_shift = 4.0, plot_theme = plot_theme)
  p_5 <- .create_density_and_qq_plot(location_shift = 5.0, plot_theme = plot_theme)

  p <- p_1$density + p_2$density + p_3$density + p_4$density + p_5$density +
    p_1$qq + p_2$qq + p_3$qq + p_4$qq + p_5$qq +
    patchwork::plot_layout(ncol = 5, heights = c(1, 0.5))

  return(p)
}


.plot_top_gear <- function(plot_theme) {

  # Plot MPG. Show Q-Q plots for:
  # - conventional transformation
  # - Raymaekers
  # - Shift-sensitive
  # - Shift-sensitive robust
  #
  # Compute empirical test p-value.


}
