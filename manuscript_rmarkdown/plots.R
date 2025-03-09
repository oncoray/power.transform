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

  # Prevent warnings due to non-standard evaluation.
  data_type <- method <- NULL

  stroke <- 0.2
  size <- 1.0

  # Process / read data.
  data <- .get_data_problematic_transformations(
    manuscript_dir = manuscript_dir
  )

  # Get annotation settings.
  annotation_settings <- get_annotation_settings(plot_theme)

  p <- ggplot2::ggplot(
    mapping = ggplot2::aes(
      y = .data$lambda,
      colour = .data$method
    )
  )
  p <- p + plot_theme
  p <- p + ggplot2::geom_hline(
    yintercept = 1.0,
    linetype = "longdash",
    colour = "grey40"
  )
  p <- p + paletteer::scale_color_paletteer_d(
    palette = "ggthemes::Tableau_10",
    name = "transformation",
    drop = FALSE
  )

  # Branch into shift and scale-specific plots
  p_shift <- p + ggplot2::scale_x_continuous(
    name = latex2exp::TeX("$\\mu$"),
    labels = scales::math_format()
  )
  p_shift <- p_shift + ggplot2::annotate(
    geom = "text",
    x = 6,
    y = 1,
    label = "expected",
    colour = annotation_settings$colour,
    family = annotation_settings$family,
    fontface = annotation_settings$face,
    size = annotation_settings$geom_text_size,
    vjust = -0.5,
    hjust = 1.0
  )

  p_scale <- p + ggplot2::scale_x_continuous(
    name = latex2exp::TeX("$\\sigma$"),
    labels = scales::math_format()
  )
  p_scale <- p_scale + ggplot2::annotate(
    geom = "text",
    x = -6,
    y = 1,
    label = "expected",
    colour = annotation_settings$colour,
    family = annotation_settings$family,
    fontface = annotation_settings$face,
    size = annotation_settings$geom_text_size,
    vjust = -0.5,
    hjust = 0.0
  )

  p_outlier <- p + ggplot2::scale_x_continuous(
    name = "d",
    labels = scales::math_format()
  )
  p_outlier <- p_outlier + ggplot2::annotate(
    geom = "text",
    x = 6,
    y = 1,
    label = "expected",
    colour = annotation_settings$colour,
    family = annotation_settings$family,
    fontface = annotation_settings$face,
    size = annotation_settings$geom_text_size,
    vjust = -0.5,
    hjust = 1.0
  )

  # Box-Cox transformation -----------------------------------------------------
  p_bc_shift <- p_shift + ggplot2::geom_point(
    data = data[method == "Box-Cox" & data_type == "shift"],
    mapping = ggplot2::aes(x = .data$y),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_bc_shift <- p_bc_shift + ggplot2::labs(title = "shift")
  p_bc_shift <- p_bc_shift + ggplot2::scale_y_continuous(
    name = latex2exp::TeX("$\\lambda$"),
    limits = c(0, 6)
  )
  p_bc_shift <- p_bc_shift + ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  p_bc_scale <- p_scale + ggplot2::geom_point(
    data = data[method == "Box-Cox" & data_type == "scale"],
    mapping = ggplot2::aes(x = .data$y),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_bc_scale <- p_bc_scale + ggplot2::labs(title = "scale")
  p_bc_scale <- p_bc_scale + ggplot2::ylim(c(0, 6))
  p_bc_scale <- p_bc_scale + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  p_bc_outlier <- p_outlier + ggplot2::geom_point(
    data = data[method == "Box-Cox" & data_type == "outlier"],
    mapping = ggplot2::aes(x = .data$y),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_bc_outlier <- p_bc_outlier + ggplot2::labs(title = "outlier")
  p_bc_outlier <- p_bc_outlier + ggplot2::ylim(c(-0.2, 1.2))
  p_bc_outlier <- p_bc_outlier + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  # Yeo-Johnson transformation -------------------------------------------------
  p_yj_shift <- p_shift + ggplot2::geom_point(
    data = data[method == "Yeo-Johnson" & data_type == "shift"],
    mapping = ggplot2::aes(x = .data$y),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_yj_shift <- p_yj_shift + ggplot2::scale_y_continuous(
    name = latex2exp::TeX("$\\lambda$"),
    limits = c(0, 6)
  )

  p_yj_scale <- p_scale + ggplot2::geom_point(
    data = data[method == "Yeo-Johnson" & data_type == "scale"],
    mapping = ggplot2::aes(x = .data$y),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_yj_scale <- p_yj_scale + ggplot2::ylim(c(0, 6))
  p_yj_scale <- p_yj_scale + ggplot2::theme(
    axis.title.y = ggplot2::element_blank()
  )

  p_yj_outlier <- p_outlier + ggplot2::geom_point(
    data = data[method == "Yeo-Johnson" & data_type == "outlier"],
    mapping = ggplot2::aes(x = .data$y),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_yj_outlier <- p_yj_outlier + ggplot2::ylim(c(-0.2, 1.2))
  p_yj_outlier <- p_yj_outlier + ggplot2::theme(
    axis.title.y = ggplot2::element_blank()
  )

  # Patch all the plots together.
  p <- p_bc_shift + p_bc_scale + p_bc_outlier +
    p_yj_shift + p_yj_scale + p_yj_outlier +
    patchwork::plot_layout(
      ncol = 3,
      guides = "collect"
    )

  return(p)
}



.plot_sensitivity_to_outlier <- function(plot_theme, manuscript_dir) {
  # Creates the plot that shows how transformation parameter lambda is affected
  # by outliers in otherwise normally distributed data..

  require(patchwork)
  require(ggplot2)

  # Prevent warnings due to non-standard evaluation.
  method <- NULL

  # Process data.
  data <- .get_shifted_outlier_plot_data(manuscript_dir = manuscript_dir)

  annotation_settings <- get_annotation_settings(plot_theme)

  p <- ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = .data$d,
      y = .data$lambda,
      colour = .data$method
    )
  )
  p <- p + plot_theme
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
    name = latex2exp::TeX("$d$"),
    labels = scales::math_format()
  )
  p <- p + ggplot2::scale_y_continuous(
    name = latex2exp::TeX("$\\lambda$"),
    limits = c(-0.2, 1.2)
  )
  p <- p + paletteer::scale_color_paletteer_d(
    name = "method",
    palette = "ggthemes::Tableau_10",
    drop = FALSE
  )

  # Box-Cox transformation -----------------------------------------------------
  p_bc <- p + ggplot2::geom_point(data = data[method == "Box-Cox"])

  # Yeo-Johnson transformation -------------------------------------------------
  p_yj <- p + ggplot2::geom_point(data = data[method == "Yeo-Johnson"])
  p_yj <- p_yj + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

  # Patch all the plots together.
  p <- p_bc + p_yj + patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  )

  return(p)
}



.plot_shifted_distributions <- function(plot_theme, manuscript_dir) {
  # Creates the plot for lambda values under location and scale shifts for MLE

  require(patchwork)
  require(ggplot2)

  # Prevent warnings due to non-standard evaluation.
  estimation_method <- distribution <- method <- version <- transformation <- NULL
  data_type <- NULL

  stroke <- 0.2
  size <- 1.0

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
  data <- .get_shifted_scaled_distribution_data(
    manuscript_dir = manuscript_dir,
    main_manuscript = TRUE
  )
  data <- data[estimation_method == "MLE"]
  data[, transformation := paste0(method, " ", version)]
  data$transformation <- factor(
    data$transformation,
    levels = c("Box-Cox conventional", "Box-Cox invariant", "Yeo-Johnson conventional", "Yeo-Johnson invariant"),
    labels = c("conventional BC", "invariant BC", "conventional YJ", "invariant YJ")
  )

  # Density plots
  p_dens <- ggplot2::ggplot(
    mapping = ggplot2::aes(x = .data$x)
  )
  p_dens <- p_dens + plot_theme
  p_dens <- p_dens + ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )

  # Lambda plots
  p <- ggplot2::ggplot(
    mapping = ggplot2::aes(
      y = .data$lambda,
      colour = .data$transformation,
      shape = .data$transformation,
    )
  )
  p <- p + plot_theme
  p <- p + ggplot2::scale_colour_discrete(
    name = "transformation",
    type = c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D"),
    drop = FALSE
  )
  p <- p + ggplot2::scale_shape_manual(
    name = "transformation",
    values = c(16, 4, 16, 4),
    drop = FALSE
  )

  # Normal distribution --------------------------------------------------------

  p_dens_normal <- p_dens + ggplot2::geom_density(data = data.table::data.table(x = x_normal))
  p_dens_normal <- p_dens_normal + ggplot2::xlim(c(-3.0, 3.0))
  p_dens_normal <- p_dens_normal + ggplot2::labs(title = "normal distribution")

  p_normal <- p + ggplot2::scale_y_continuous(
    name = latex2exp::TeX("$\\lambda$"),
    limits = c(0.0, 6.0)
  )

  p_normal_shift <- p_normal + ggplot2::scale_x_continuous(
    name = latex2exp::TeX("$d_{shift}$"),
    labels = scales::math_format()
  )

  p_normal_scale <- p_normal + ggplot2::scale_x_continuous(
    name = latex2exp::TeX("$d_{scale}$"),
    labels = scales::math_format()
  )

  ## Box-Cox transformation ----------------------------------------------------
  p_normal_shift_bc <- p_normal_shift + ggplot2::geom_point(
    data = data[distribution == "normal" & method == "Box-Cox" & data_type == "shifted"],
    mapping = ggplot2::aes(x = .data$shift),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_normal_shift_bc <- p_normal_shift_bc + ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  p_normal_scale_bc <- p_normal_scale + ggplot2::geom_point(
    data = data[distribution == "normal" & method == "Box-Cox" & data_type == "scaled"],
    mapping = ggplot2::aes(x = .data$scale),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_normal_scale_bc <- p_normal_scale_bc + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  ## Yeo-Johnson transformation ------------------------------------------------

  p_normal_shift_yj <- p_normal_shift + ggplot2::geom_point(
    data = data[distribution == "normal" & method == "Yeo-Johnson" & data_type == "shifted"],
    mapping = ggplot2::aes(x = .data$shift),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )

  p_normal_scale_yj <- p_normal_scale + ggplot2::geom_point(
    data = data[distribution == "normal" & method == "Yeo-Johnson" & data_type == "scaled"],
    mapping = ggplot2::aes(x = .data$scale),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_normal_scale_yj <- p_normal_scale_yj + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

  # Right-skewed distribution --------------------------------------------------

  p_dens_right <- p_dens + ggplot2::geom_density(data = data.table::data.table(x = x_right_skewed))
  p_dens_right <- p_dens_right + ggplot2::xlim(c(-3.0, 10.0))
  p_dens_right <- p_dens_right + ggplot2::labs(title = "right-skewed distribution")

  p_right <- p + ggplot2::scale_y_continuous(
    name = latex2exp::TeX("$\\lambda$"),
    limits = c(-4.5, 1.5)
  )

  p_right_shift <- p_right + ggplot2::scale_x_continuous(
    name = latex2exp::TeX("$d_{shift}$"),
    labels = scales::math_format()
  )

  p_right_scale <- p_right + ggplot2::scale_x_continuous(
    name = latex2exp::TeX("$d_{scale}$"),
    labels = scales::math_format()
  )

  ## Box-Cox transformation ----------------------------------------------------
  p_right_shift_bc <- p_right_shift + ggplot2::geom_point(
    data = data[distribution =="right-skewed" & method == "Box-Cox" & data_type == "shifted"],
    mapping = ggplot2::aes(x = .data$shift),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_right_shift_bc <- p_right_shift_bc + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  p_right_scale_bc <- p_right_scale + ggplot2::geom_point(
    data = data[distribution == "right-skewed" & method == "Box-Cox" & data_type == "scaled"],
    mapping = ggplot2::aes(x = .data$scale),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_right_scale_bc <- p_right_scale_bc + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  ## Yeo-Johnson transformation ------------------------------------------------

  p_right_shift_yj <- p_right_shift + ggplot2::geom_point(
    data = data[distribution == "right-skewed" & method == "Yeo-Johnson" & data_type == "shifted"],
    mapping = ggplot2::aes(x = .data$shift),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_right_shift_yj <- p_right_shift_yj + ggplot2::theme(
    axis.title.y = ggplot2::element_blank()
  )

  p_right_scale_yj <- p_right_scale + ggplot2::geom_point(
    data = data[distribution == "right-skewed" & method == "Yeo-Johnson" & data_type == "scaled"],
    mapping = ggplot2::aes(x = .data$scale),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_right_scale_yj <- p_right_scale_yj + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )


  # Left-skewed distribution --------------------------------------------------

  p_dens_left <- p_dens + ggplot2::geom_density(data = data.table::data.table(x = x_left_skewed))
  p_dens_left <- p_dens_left + ggplot2::xlim(c(-10.0, 3.0))
  p_dens_left <- p_dens_left + ggplot2::labs(title = "left-skewed distribution")

  p_left <- p + ggplot2::scale_y_continuous(
    name = latex2exp::TeX("$\\lambda$"),
    limits = c(0.0, 6.0)
  )

  p_left_shift <- p_left + ggplot2::scale_x_continuous(
    name = latex2exp::TeX("$d_{shift}$"),
    labels = scales::math_format()
  )

  p_left_scale <- p_left + ggplot2::scale_x_continuous(
    name = latex2exp::TeX("$d_{scale}$"),
    labels = scales::math_format()
  )

  ## Box-Cox transformation ----------------------------------------------------
  p_left_shift_bc <- p_left_shift + ggplot2::geom_point(
    data = data[distribution =="left-skewed" & method == "Box-Cox" & data_type == "shifted"],
    mapping = ggplot2::aes(x = .data$shift),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_left_shift_bc <- p_left_shift_bc + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  p_left_scale_bc <- p_left_scale + ggplot2::geom_point(
    data = data[distribution == "left-skewed" & method == "Box-Cox" & data_type == "scaled"],
    mapping = ggplot2::aes(x = .data$scale),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_left_scale_bc <- p_left_scale_bc + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  ## Yeo-Johnson transformation ------------------------------------------------

  p_left_shift_yj <- p_left_shift + ggplot2::geom_point(
    data = data[distribution == "left-skewed" & method == "Yeo-Johnson" & data_type == "shifted"],
    mapping = ggplot2::aes(x = .data$shift),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_left_shift_yj <- p_left_shift_yj + ggplot2::theme(
    axis.title.y = ggplot2::element_blank()
  )

  p_left_scale_yj <- p_left_scale + ggplot2::geom_point(
    data = data[distribution == "left-skewed" & method == "Yeo-Johnson" & data_type == "scaled"],
    mapping = ggplot2::aes(x = .data$scale),
    stroke = stroke,
    size = size,
    show.legend = TRUE
  )
  p_left_scale_yj <- p_left_scale_yj + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )


  # Patch all the plots together.
  p <-
    (p_dens_normal | p_dens_right | p_dens_left ) /
    (p_normal_shift_bc | p_normal_scale_bc | p_right_shift_bc | p_right_scale_bc | p_left_shift_bc | p_left_scale_bc) /
    (p_normal_shift_yj | p_normal_scale_yj | p_right_shift_yj | p_right_scale_yj | p_left_shift_yj | p_left_scale_yj) +
    patchwork::plot_layout(
      heights = c(0.5, 1.0, 1.0),
      guides = "collect",
    )

  return(p)
}



.plot_shifted_distributions_appendix <- function(plot_theme, manuscript_dir) {
  # Creates the plot for lambda values under location shift for all optimisation
  # criteria.

  # Prevent warnings due to non-standard evaluation.
  distribution <- method <- data_type <- version <- NULL

  # Update plot theme
  plot_theme$plot.tag.position <- "top"
  plot_theme$plot.tag$margin <- grid::unit(c(0, 0, 5, 0), "points")
  plot_theme$plot.tag$face <- "bold"

  ..create_plot <- function(
    p,
    distr,
    meth,
    ver,
    type,
    first_col = FALSE,
    last_row = FALSE,
    stroke = 0.2,
    size = 1.0
  ) {

    if (type == "shifted") {
      p <- p + ggplot2::geom_point(
        data = data[
          distribution == distr &
            method == meth &
            version == ver &
            data_type == type
        ],
        mapping = ggplot2::aes(x = .data$shift),
        stroke = stroke,
        size = size
      )

      p <- p + ggplot2::scale_x_continuous(
        name = latex2exp::TeX("$d_{shift}$"),
        labels = scales::math_format()
      )

    } else if (type == "scaled") {
      p <- p + ggplot2::geom_point(
        data = data[
          distribution == distr &
            method == meth &
            version == ver &
            data_type == type
        ],
        mapping = ggplot2::aes(x = .data$scale),
        stroke = stroke,
        size = size
      )

      p <- p + ggplot2::scale_x_continuous(
        name = latex2exp::TeX("$d_{scale}$"),
        labels = scales::math_format()
      )
    }

    if (type == "scaled") {
      p <- p + ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
    }
    if (!first_col) {
      p <- p + ggplot2::theme(
        axis.title.y = ggplot2::element_blank()
      )
    }
    if (!last_row) {
      p <- p + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    }

    return(p)
  }

  # Process / read data.
  data <- .get_shifted_scaled_distribution_data(
    manuscript_dir = manuscript_dir,
    main_manuscript = FALSE
  )

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

  # Density plots
  p_dens <- ggplot2::ggplot(
    mapping = ggplot2::aes(x = .data$x)
  )
  p_dens <- p_dens + plot_theme
  p_dens <- p_dens + ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )

  # Lambda plots
  p <- ggplot2::ggplot(
    mapping = ggplot2::aes(
      y = .data$lambda,
      colour = .data$estimation_method
    )
  )
  p <- p + plot_theme
  p <- p + paletteer::scale_color_paletteer_d(
    name = "optimisation criterion",
    palette = "ggthemes::Tableau_10",
    drop = FALSE
  )

  # Normal distribution --------------------------------------------------------

  p_dens_normal <- p_dens + ggplot2::geom_density(data = data.table::data.table(x = x_normal))
  p_dens_normal <- p_dens_normal + ggplot2::xlim(c(-3.0, 3.0))
  p_dens_normal <- p_dens_normal + ggplot2::labs(title = "normal distribution")

  p_normal <- p + ggplot2::scale_y_continuous(
    name = latex2exp::TeX("$\\lambda$"),
    limits = c(-4.0, 6.0)
  )

  ## Box-Cox transformation ----------------------------------------------------
  p_normal_shift_conv_bc <- ..create_plot(p_normal, "normal", "Box-Cox", "conventional", "shifted", first_col = TRUE)
  p_normal_shift_invar_bc <- ..create_plot(p_normal, "normal", "Box-Cox", "invariant", "shifted", first_col = TRUE)
  p_normal_scale_conv_bc <- ..create_plot(p_normal, "normal", "Box-Cox", "conventional", "scaled")
  p_normal_scale_invar_bc <- ..create_plot(p_normal, "normal", "Box-Cox", "invariant", "scaled")

  ## Yeo-Johnson transformation ------------------------------------------------
  p_normal_shift_conv_yj <- ..create_plot(p_normal, "normal", "Yeo-Johnson", "conventional", "shifted", first_col = TRUE)
  p_normal_shift_invar_yj <- ..create_plot(p_normal, "normal", "Yeo-Johnson", "invariant", "shifted", first_col = TRUE, last_row = TRUE)
  p_normal_scale_conv_yj <- ..create_plot(p_normal, "normal", "Yeo-Johnson", "conventional", "scaled")
  p_normal_scale_invar_yj <- ..create_plot(p_normal, "normal", "Yeo-Johnson", "invariant", "scaled", last_row = TRUE)


  # Right-skewed distribution --------------------------------------------------

  p_dens_right <- p_dens + ggplot2::geom_density(data = data.table::data.table(x = x_right_skewed))
  p_dens_right <- p_dens_right + ggplot2::xlim(c(-3.0, 10.0))
  p_dens_right <- p_dens_right + ggplot2::labs(title = "right-skewed distribution")

  p_right <- p + ggplot2::scale_y_continuous(
    name = latex2exp::TeX("$\\lambda$"),
    limits = c(-4.5, 1.5)
  )

  ## Box-Cox transformation ----------------------------------------------------
  p_right_shift_conv_bc <- ..create_plot(p_right, "right-skewed", "Box-Cox", "conventional", "shifted")
  p_right_shift_invar_bc <- ..create_plot(p_right, "right-skewed", "Box-Cox", "invariant", "shifted")
  p_right_scale_conv_bc <- ..create_plot(p_right, "right-skewed", "Box-Cox", "conventional", "scaled")
  p_right_scale_invar_bc <- ..create_plot(p_right, "right-skewed", "Box-Cox", "invariant", "scaled")

  ## Yeo-Johnson transformation ------------------------------------------------
  p_right_shift_conv_yj <- ..create_plot(p_right, "right-skewed", "Yeo-Johnson", "conventional", "shifted")
  p_right_shift_invar_yj <- ..create_plot(p_right, "right-skewed", "Yeo-Johnson", "invariant", "shifted", last_row = TRUE)
  p_right_scale_conv_yj <- ..create_plot(p_right, "right-skewed", "Yeo-Johnson", "conventional", "scaled")
  p_right_scale_invar_yj <- ..create_plot(p_right, "right-skewed", "Yeo-Johnson", "invariant", "scaled", last_row = TRUE)


  # Left-skewed distribution --------------------------------------------------

  p_dens_left <- p_dens + ggplot2::geom_density(data = data.table::data.table(x = x_left_skewed))
  p_dens_left <- p_dens_left + ggplot2::xlim(c(-10.0, 3.0))
  p_dens_left <- p_dens_left + ggplot2::labs(title = "left-skewed distribution")

  p_left <- p + ggplot2::scale_y_continuous(
    name = latex2exp::TeX("$\\lambda$"),
    limits = c(0.0, 6.0)
  )

  ## Box-Cox transformation ----------------------------------------------------
  p_left_shift_conv_bc <- ..create_plot(p_left, "left-skewed", "Box-Cox", "conventional", "shifted")
  p_left_shift_invar_bc <- ..create_plot(p_left, "left-skewed", "Box-Cox", "invariant", "shifted")
  p_left_scale_conv_bc <- ..create_plot(p_left, "left-skewed", "Box-Cox", "conventional", "scaled")
  p_left_scale_invar_bc <- ..create_plot(p_left, "left-skewed", "Box-Cox", "invariant", "scaled")

  ## Yeo-Johnson transformation ------------------------------------------------
  p_left_shift_conv_yj <- ..create_plot(p_left, "left-skewed", "Yeo-Johnson", "conventional", "shifted")
  p_left_shift_invar_yj <- ..create_plot(p_left, "left-skewed", "Yeo-Johnson", "invariant", "shifted", last_row = TRUE)
  p_left_scale_conv_yj <- ..create_plot(p_left, "left-skewed", "Yeo-Johnson", "conventional", "scaled")
  p_left_scale_invar_yj <- ..create_plot(p_left, "left-skewed", "Yeo-Johnson", "invariant", "scaled", last_row = TRUE)

  # Add tags
  p_normal_shift_conv_bc <- p_normal_shift_conv_bc + ggplot2::labs(tag = "conventional Box-Cox")
  p_normal_shift_invar_bc <- p_normal_shift_invar_bc + ggplot2::labs(tag = "invariant Box-Cox")
  p_normal_shift_conv_yj <- p_normal_shift_conv_yj + ggplot2::labs(tag = "conventional Yeo-Johnson")
  p_normal_shift_invar_yj <- p_normal_shift_invar_yj + ggplot2::labs(tag = "invariant Yeo-Johnson")

  # Patch all the plots together.
  p <-
    (p_dens_normal | p_dens_right | p_dens_left ) /
    (p_normal_shift_conv_bc  | p_normal_scale_conv_bc  | p_right_shift_conv_bc  | p_right_scale_conv_bc  | p_left_shift_conv_bc  | p_left_scale_conv_bc) /
    (p_normal_shift_invar_bc | p_normal_scale_invar_bc | p_right_shift_invar_bc | p_right_scale_invar_bc | p_left_shift_invar_bc | p_left_scale_invar_bc) /
    (p_normal_shift_conv_yj  | p_normal_scale_conv_yj  | p_right_shift_conv_yj  | p_right_scale_conv_yj  | p_left_shift_conv_yj  | p_left_scale_conv_yj) /
    (p_normal_shift_invar_yj | p_normal_scale_invar_yj | p_right_shift_invar_yj | p_right_scale_invar_yj | p_left_shift_invar_yj | p_left_scale_invar_yj) +
    patchwork::plot_layout(
      heights = c(0.5, 1.0, 1.0, 1.0, 1.0),
      guides = "collect",
    )

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
  p_step <- p_step + ggplot2::labs(title = "step")

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
  p_triangle <- p_triangle + ggplot2::labs(title = "triangle")
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
  p_cosine <- p_cosine + ggplot2::labs(title = "tapered cosine")
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
    drop = FALSE,
    name = "transformation"
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
  p <- p + ggplot2::scale_y_sqrt(name = latex2exp::TeX("$| \\lambda_r - \\lambda_0 |$"))
  p <- p + ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
  p <- p + ggplot2::xlab("weighting method")

  return(p)
}



.plot_residuals <- function(plot_theme, manuscript_dir) {
  # Prevent warnings due to non-standard evaluation.
  has_outliers <- residual <- threshold <- residual_error <- method <- p <- NULL

  if (!file.exists(file.path(manuscript_dir, "manuscript_residual_error.RDS"))) {
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
      file.path(manuscript_dir, "manuscript_residual_error.RDS")
    )
  } else {
    data <- readRDS(file.path(manuscript_dir, "manuscript_residual_error.RDS"))

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
    name = "percentile\n(Box-Cox)",
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
    name = "percentile\n(Yeo-Johnson)",
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
  p_bc <- p_bc + ggplot2::xlab(latex2exp::TeX("test statistic $\\tau_{ecn}$"))
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
  p_yj <- p_yj + ggplot2::xlab(latex2exp::TeX("test statistic $\\tau_{ecn}$"))
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
  outlier_id <- test_statistic <- NULL

  file_name <- file.path(manuscript_dir, "type_1_error_rate_plot_appendix.RDS")

  if (!file.exists(file_name)) {
    # Get data. We need three datasets:
    # - Dataset presented in main manuscript.
    # -
    data <- .get_goodness_of_fit_data(manuscript_dir = manuscript_dir)
    data_no_outlier <- data.table::copy(data[outlier_id == 1L])
    data_no_outlier$method <- factor(
      data_no_outlier$method,
      levels = c("Box-Cox", "Yeo-Johnson"),
      labels = c("Box-Cox (no outlier)", "Yeo-Johnson (no outlier)")
    )
    data_normal <- .get_goodness_of_fit_data_normal_only(manuscript_dir = manuscript_dir)
    data <- data.table::rbindlist(list(data, data_no_outlier, data_normal))


    # Process data
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
          "mare" = c(0.0, 0.0, 0.0, 0.0, 0.0),
          "method" = c("Box-Cox", "Yeo-Johnson", "Box-Cox (no outlier)", "Yeo-Johnson (no outlier)", "none"),
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

    saveRDS(data, file_name)

  } else {
    data <- readRDS(file_name)
  }

  data$method <- factor(
    data$method,
    c("Box-Cox", "Box-Cox (no outlier)", "Yeo-Johnson", "Yeo-Johnson (no outlier)", "none")
  )

  p <- ggplot2::ggplot(
    data = data[kappa == "80 %", ],
    mapping = ggplot2::aes(
      x = mare,
      y = rejected,
      colour = method
    )
  )
  p <- p + plot_theme
  p <- p + ggplot2::geom_step()
  p <- p + ggplot2::scale_colour_discrete(
    name = "type",
    type = c(
      "Box-Cox" = "#4E79A7",
      "Box-Cox (no outlier)" = "#A0CBE8",
      "Yeo-Johnson" = "#F28E2B",
      "Yeo-Johnson (no outlier)" = "#FFBE7D",
      "none" = "#59A14F"
    )
  )
  p <- p + ggplot2::xlab(latex2exp::TeX("test statistic $\\tau_{ecn}$"))
  p <- p + ggplot2::ylab("type I error rate")
  p <- p + ggplot2::coord_cartesian(xlim = c(0, 0.25))

  # Process test data.
  significance_values <- c(0.50, 0.20, 0.10, 0.05, 0.02, 0.01, 0.001)

  test_data <- data[
    kappa == "80 %",
    list(
      "test_statistic" = stats::spline(
        x = rejected,
        y = mare,
        method = "hyman",
        xout = significance_values
      )$y,
      "alpha" = significance_values
    ),
    by = "method"
  ]

  test_data$alpha <- factor(test_data$alpha, levels=significance_values)
  test_data[, "test_statistic" := round(test_statistic, digits=3)]

  test_data <- data.table::dcast(
    data = test_data,
    method ~ alpha,
    value.var = "test_statistic"
  )

  return(list("p" = p, "data" = test_data))
}

#
#
#
# .plot_type_1_error_rate_appendix <- function(plot_theme, manuscript_dir) {
#   # Prevent warnings due to non-standard evaluation.
#   residual <- rejected <- p <- n <- mare <- method <- NULL
#
#   # Get data
#   data <- .get_test_statistics_data_appendix(manuscript_dir = manuscript_dir)
#
#   central_width <- c(0.60, 0.70, 0.80, 0.90, 0.95, 1.00)
#
#   data$kappa <- factor(
#     x = data$kappa,
#     levels = central_width,
#     labels = c("60 %", "70 %", "80 %", "90 %", "95 %", "100 %")
#   )
#
#   data <- data[kappa == "80 %"]
#
#   plot_start_list <- list()
#   ii <- 1L
#   for (method in levels(data$method)) {
#     for (method_tag in levels(data$method_tag)) {
#       for (kappa in levels(data$kappa)) {
#         plot_start_list[[ii]] <- data.table::data.table(
#           mare = 0.0,
#           method = method,
#           method_tag = method_tag,
#           n = 0L,
#           rejected = 1.0,
#           kappa = kappa
#         )
#
#         ii <- ii + 1L
#       }
#     }
#   }
#
#   data <- data.table::rbindlist(c(list(data), plot_start_list))
#
#   p_bc <- ggplot2::ggplot(
#     data = data[method == "Box-Cox"],
#     mapping = ggplot2::aes(
#       x = mare,
#       y = rejected,
#       colour = method_tag
#     )
#   )
#   p_bc <- p_bc + plot_theme
#   p_bc <- p_bc + ggplot2::geom_step()
#   p_bc <- p_bc + ggplot2::scale_colour_discrete(
#     name = "power transform variant\n(Box-Cox)",
#     type = c(
#       "robust, shift sens. (MLE)" = "#bacbde",
#       "conventional (MLE)" = "#98b2cd",
#       "robust, shift sens. (C-vM)" = "#42648a",
#       "conventional (C-vM)" = "#324b67"
#     )
#   )
#   p_bc <- p_bc + ggplot2::xlab("test statistic")
#   p_bc <- p_bc + ggplot2::ylab("type I error rate")
#   p_bc <- p_bc + ggplot2::coord_cartesian(xlim = c(0, 0.15))
#
#   p_yj <- ggplot2::ggplot(
#     data = data[method == "Yeo-Johnson"],
#     mapping = ggplot2::aes(
#       x = mare,
#       y = rejected,
#       colour = method_tag
#     )
#   )
#   p_yj <- p_yj + plot_theme
#   p_yj <- p_yj + ggplot2::geom_step()
#   p_yj <- p_yj + ggplot2::scale_colour_discrete(
#     name = "power transform variant\n(Yeo-Johnson)",
#     type = c(
#       "robust, shift sens. (MLE)" = "#f9c59f",
#       "conventional (MLE)" = "#f6a76f",
#       "robust, shift sens. (C-vM)" = "#c0570c",
#       "conventional (C-vM)" = "#904109"
#     )
#   )
#   p_yj <- p_yj + ggplot2::xlab("test statistic")
#   p_yj <- p_yj + ggplot2::coord_cartesian(xlim = c(0, 0.15))
#   p_yj <- p_yj + ggplot2::theme(
#     axis.title.y = ggplot2::element_blank(),
#     axis.text.y = ggplot2::element_blank()
#   )
#
#   p <- p_bc + p_yj + patchwork::plot_layout(ncol = 2, guides = "collect")
#
#   return(p)
# }


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
      vjust = -0.5,
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

    p_qq <- p_qq + ggplot2::xlab("expected quantile")
    p_qq <- p_qq + ggplot2::ylab("observed quantile")

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
    patchwork::plot_layout(ncol = 5, heights = c(1.3, 1))

  return(p)
}


.plot_lung_cancer_age <- function(plot_theme, lambda_limit = NULL) {
  require(survival)

  return(
    ..plot_real_data(
      plot_theme = plot_theme,
      x = survival::lung$age,
      lambda_limit = NULL
    )
  )
}



.plot_penguin_body_mass <- function(plot_theme, lambda_limit = NULL) {
  require(modeldata)

  x <- modeldata::penguins$body_mass_g
  x <- x[is.finite(x)]

  return(
    ..plot_real_data(
      plot_theme = plot_theme,
      x = x,
      lambda_limit = NULL
    )
  )
}



.plot_ames_latitude <- function(plot_theme, lambda_limit = NULL) {
  require(modeldata)

  return(
    ..plot_real_data(
      plot_theme = plot_theme,
      x = modeldata::ames$Latitude,
      lambda_limit = NULL
    )
  )
}



.plot_top_gear <- function(plot_theme) {
  require(robustHD)

  # Load TopGear data into the current environment
  utils::data("TopGear", package = "robustHD", envir = environment())
  data <- data.table::as.data.table(TopGear)

  # Extract MPG data.
  x <- data$MPG
  x <- x[!is.na(x)]

  return(
    ..plot_real_data(
      plot_theme = plot_theme,
      x = x
    )
  )
}



.plot_ischemic_stroke <- function(plot_theme) {
  require(modeldata)

  x <- modeldata::ischemic_stroke$max_max_wall_thickness

  return(
    ..plot_real_data(
      plot_theme = plot_theme,
      x = x
    )
  )
}



..plot_real_data <- function(plot_theme, x, lambda_limit = c(-4.0, 6.0)) {
  # Show density and Q-Q plots for:
  # - Original data
  # - conventional transformation
  # - Raymaekers
  # - Shift-sensitive
  # - Shift-sensitive robust
  #
  # Compute empirical test p-value.

  # Prevent undue warnings.
  transformation <- z_observed <- NULL

  plot_limits <- c(-3.0, 3.0)

  ..parse_data <- function(
    x,
    transformer,
    name
  ) {
    transformed_data <- power.transform::power_transform(
      x = x,
      transformer = transformer
    )

    summary_data <- data.table::data.table(
      "mu" = mean(transformed_data),
      "sigma" = sd(transformed_data),
      "lambda" = power.transform::get_lambda(transformer),
      "shift" = power.transform::get_shift(transformer),
      "scale" = power.transform::get_scale(transformer),
      "transformation" = name
    )

    transformed_data <- data.table::data.table(
      x = transformed_data,
      x_standardised = (transformed_data - stats::median(transformed_data)) / stats::IQR(transformed_data),
      transformation = name
    )

    residual_data <- power.transform::get_residuals(
      x = x,
      transformer = transformer
    )

    residual_data[, "transformation" := name]

    central_normality_data <- data.table::data.table(
      p_value = power.transform::assess_transformation(
        x = x,
        transformer = transformer,
        verbose = FALSE
      ),
      transformation = name
    )

    return(list(
      "transformed_data" = transformed_data,
      "residual_data" = residual_data,
      "summary_data" = summary_data,
      "central_normality_data" = central_normality_data
    ))
  }

  # No transformer
  no_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = "none"
  )

  # Conventional transformer
  yj_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = "yeo_johnson",
    robust = FALSE,
    invariant = FALSE,
    lambda = lambda_limit
  )

  # Raymaekers robust transformer
  yj_rr_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = "yeo_johnson",
    robust = TRUE,
    invariant = FALSE,
    estimation_method = "raymaekers_robust",
    lambda = lambda_limit
  )

  # Shift-sensitive transformer
  yj_shift_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = "yeo_johnson",
    robust = FALSE,
    invariant = TRUE,
    lambda = lambda_limit
  )

  # Robust shift-sensitive transformer
  yj_shift_robust_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = "yeo_johnson",
    robust = TRUE,
    invariant = TRUE,
    lambda = lambda_limit
  )

  transformer_labels <- c(
    "none", "conventional", "Raymaekers-Rousseeuw", "invariant", "robust invariant"
  )

  data <- mapply(
    FUN = ..parse_data,
    transformer = list(
      no_transformer,
      yj_transformer,
      yj_rr_transformer,
      yj_shift_transformer,
      yj_shift_robust_transformer
    ),
    name = transformer_labels,
    MoreArgs = list("x" = x),
    SIMPLIFY = FALSE
  )
  names(data) <- transformer_labels

  # Parse transformed data
  transformed_data <- data.table::rbindlist(
    lapply(data, function(x) (x$transformed_data)),
    use.names = TRUE
  )
  transformed_data$transformation <- factor(
    transformed_data$transformation,
    levels = transformer_labels
  )

  # Parse residual data
  residual_data <- data.table::rbindlist(
    lapply(data, function(x) (x$residual_data)),
    use.names = TRUE
  )
  residual_data$transformation <- factor(
    residual_data$transformation,
    levels = transformer_labels
  )
  residual_data[, ":="("z_observed_truncated" = z_observed, "truncated" = FALSE)]
  residual_data[z_observed < plot_limits[1], ":="("z_observed_truncated" = plot_limits[1], "truncated" = TRUE)]
  residual_data[z_observed > plot_limits[2], ":="("z_observed_truncated" = plot_limits[2], "truncated" = TRUE)]

  # Set robust tag
  residual_data[]

  # Parse central normality data.
  central_normality_data <- data.table::rbindlist(
    lapply(data, function(x) (x$central_normality_data)),
    use.names = TRUE
  )

  # Quantile-quantile plots.
  p_qq <- ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = .data$z_expected,
      y = .data$z_observed_truncated,
      colour = .data$transformation,
      size = .data$transformation,
      shape = .data$truncated
    ))
  p_qq <- p_qq + plot_theme
  p_qq <- p_qq + ggplot2::scale_colour_manual(
    name = "transformation",
    values = c("#111111", "#8cd17d", "#59a14f", "#ff9d9a", "#e15759"),
    breaks = transformer_labels,
    drop = FALSE
  )
  p_qq <- p_qq + ggplot2::scale_size_manual(
    name = "transformation",
    values = c(1.2, 1.2, 0.6, 1.2, 0.6),
    breaks = transformer_labels,
    drop = FALSE
  )
  p_qq <- p_qq + ggplot2::scale_shape_manual(
    values = c(16, 4),
    labels = c(FALSE, TRUE),
    guide = "none"
  )
  p_qq <- p_qq + ggplot2::coord_cartesian(
    xlim = plot_limits,
    ylim = plot_limits
  )
  p_qq <- p_qq + ggplot2::geom_abline(
    intercept = 0.0,
    slope = 1.0
  )
  p_qq <- p_qq + ggplot2::xlab("expected quantile")
  p_qq <- p_qq + ggplot2::ylab("observed quantile")

  # Density plots.
  p_d <- ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = .data$x_standardised,
      y = ggplot2::after_stat(scaled),
      colour = .data$transformation
    ))
  p_d <- p_d + plot_theme
  p_d <- p_d + ggplot2::scale_colour_manual(
    name = "transformation",
    values = c("#111111", "#8cd17d", "#59a14f", "#ff9d9a", "#e15759"),
    breaks = transformer_labels,
    drop = FALSE,
    guide = "none"
  )
  p_d <- p_d + ggplot2::xlim(c(-3, 3))
  p_d <- p_d + ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )

  # Without transformation -----------------------------------------------------
  p0_qq <- p_qq + ggplot2::geom_point(
    data = residual_data[transformation == "none"]
  )
  p0_d <- p_d + ggplot2::geom_density(
    data = transformed_data[transformation == "none"]
  )
  p0_d <- p0_d + ggplot2::labs(title = "none")

  # Standard transformation ----------------------------------------------------
  p1_qq <- p_qq + ggplot2::geom_point(
    data  = residual_data[transformation %in% c("conventional", "Raymaekers-Rousseeuw")]
  )
  p1_qq <- p1_qq + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )
  p1_d <- p_d + ggplot2::geom_density(
    data = transformed_data[transformation %in% c("conventional", "Raymaekers-Rousseeuw")]
  )
  p1_d <- p1_d + ggplot2::labs(title = "conventional")

  # Shift-sensitive transformation ---------------------------------------------
  p2_qq <- p_qq + ggplot2::geom_point(
    data  = residual_data[transformation %in% c("invariant", "robust invariant")]
  )
  p2_qq <- p2_qq + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )
  p2_d <- p_d + ggplot2::geom_density(
    data = transformed_data[transformation %in% c("invariant", "robust invariant")]
  )
  p2_d <- p2_d + ggplot2::labs(title = "invariant")

  p <- p0_d + p1_d + p2_d + p0_qq + p1_qq + p2_qq + patchwork::plot_layout(
    ncol = 3,
    guides = "collect",
    heights = c(0.3, 1.0)
  )

  return(list("data" = data, "plot"= p))
}



.plot_experimental_invariance <- function(plot_theme, lambda_limit = NULL, method = "yeo_johnson") {

  # Age
  p_age <- ..plot_combined_real_data(
    plot_theme = plot_theme,
    x = survival::lung$age,
    data_name = "age",
    lambda_limit = lambda_limit,
    method = method,
    first_set = TRUE
  )

  # Penguin body mass
  x <- modeldata::penguins$body_mass_g
  x <- x[is.finite(x)]
  p_pen <- ..plot_combined_real_data(
    plot_theme = plot_theme,
    x = x,
    data_name = "penguin body mass",
    lambda_limit = lambda_limit,
    method = method
  )

  # Housing latitude
  p_lat <- ..plot_combined_real_data(
    plot_theme = plot_theme,
    x = modeldata::ames$Latitude,
    data_name = "latitude",
    lambda_limit = lambda_limit,
    method = method,
    last_set = TRUE
  )

  p <- p_age$density[[1]] + p_age$density[[2]] + p_age$density[[3]] +
    p_age$qq[[1]] + p_age$qq[[2]] + p_age$qq[[3]] +
    p_pen$density[[1]] + p_pen$density[[2]] + p_pen$density[[3]] +
    p_pen$qq[[1]] + p_pen$qq[[2]] + p_pen$qq[[3]] +
    p_lat$density[[1]] + p_lat$density[[2]] + p_lat$density[[3]] +
    p_lat$qq[[1]] + p_lat$qq[[2]] + p_lat$qq[[3]] +
    patchwork::plot_layout(
      ncol = 3,
      guides = "collect",
      heights = c(0.2, 1.0, 0.2, 1.0, 0.2, 1.0)
    )

  return(list(
    "data" = list("age" = p_age$data, "penguin_body_mass" = p_pen$data, "latitude" = p_lat$data),
    "plot"= p
  ))
}



.plot_experimental_outlier_robustness <- function(plot_theme, lambda_limit = NULL, method = "yeo_johnson") {
  # Top gear data
  utils::data("TopGear", package = "robustHD", envir = environment())
  data <- data.table::as.data.table(TopGear)
  x <- data$MPG
  x <- x[!is.na(x)]
  p_top <- ..plot_combined_real_data(
    plot_theme = plot_theme,
    x = x,
    data_name = "fuel efficiency",
    lambda_limit = lambda_limit,
    method = method,
    first_set = TRUE
  )

  # Ischemic stroke
  p_mwt <- ..plot_combined_real_data(
    plot_theme = plot_theme,
    x = modeldata::ischemic_stroke$max_max_wall_thickness,
    data_name = "arterial wall thickness",
    lambda_limit = lambda_limit,
    method = method,
    last_set = TRUE
  )

  p <- p_top$density[[1]] + p_top$density[[2]] + p_top$density[[3]] +
    p_top$qq[[1]] + p_top$qq[[2]] + p_top$qq[[3]] +
    p_mwt$density[[1]] + p_mwt$density[[2]] + p_mwt$density[[3]] +
    p_mwt$qq[[1]] + p_mwt$qq[[2]] + p_mwt$qq[[3]] +
    patchwork::plot_layout(
      ncol = 3,
      guides = "collect",
      heights = c(0.2, 1.0, 0.2, 1.0)
    )

  return(list(
    "data" = list("fuel_efficiency" = p_top$data, "arterial_wall_thickness" = p_mwt$data),
    "plot"= p
  ))
}



..plot_combined_real_data <- function(
    plot_theme,
    x,
    data_name,
    lambda_limit = c(-4.0, 6.0),
    method = "yeo_johnson",
    first_set = FALSE,
    last_set = FALSE
) {
  # Show density and Q-Q plots for:
  # - Original data
  # - conventional transformation
  # - Raymaekers
  # - location- and shift-invariant transformation
  # - robust location- and shift-invariant transformation
  #
  # Compute empirical test p-value.

  # Prevent undue warnings.
  transformation <- z_observed <- NULL

  plot_limits <- c(-3.0, 3.0)

  ..parse_data <- function(
    x,
    transformer,
    name
  ) {
    transformed_data <- power.transform::power_transform(
      x = x,
      transformer = transformer
    )

    summary_data <- data.table::data.table(
      "mu" = mean(transformed_data),
      "sigma" = sd(transformed_data),
      "lambda" = power.transform::get_lambda(transformer),
      "shift" = power.transform::get_shift(transformer),
      "scale" = power.transform::get_scale(transformer),
      "transformation" = name
    )

    transformed_data <- data.table::data.table(
      x = transformed_data,
      x_standardised = (transformed_data - stats::median(transformed_data)) / stats::IQR(transformed_data),
      transformation = name
    )

    residual_data <- power.transform::get_residuals(
      x = x,
      transformer = transformer
    )

    residual_data[, "transformation" := name]

    central_normality_data <- data.table::data.table(
      p_value = power.transform::assess_transformation(
        x = x,
        transformer = transformer,
        verbose = FALSE
      ),
      transformation = name
    )

    return(list(
      "transformed_data" = transformed_data,
      "residual_data" = residual_data,
      "summary_data" = summary_data,
      "central_normality_data" = central_normality_data
    ))
  }

  # No transformer
  no_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = "none"
  )

  # Conventional transformer
  conv_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = method,
    robust = FALSE,
    invariant = FALSE,
    lambda = lambda_limit
  )

  # Raymaekers robust transformer
  rr_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = method,
    robust = TRUE,
    invariant = FALSE,
    estimation_method = "raymaekers_robust",
    lambda = lambda_limit
  )

  # Location- and scale-invariant transformer
  invar_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = method,
    robust = FALSE,
    invariant = TRUE,
    lambda = lambda_limit
  )

  # Robust location- and scale-invariant transformer
  invar_robust_transformer <- power.transform::find_transformation_parameters(
    x = x,
    method = method,
    robust = TRUE,
    invariant = TRUE,
    lambda = lambda_limit
  )

  transformer_labels <- c(
    "none", "conventional", "Raymaekers-Rousseeuw", "invariant", "robust invariant"
  )

  data <- mapply(
    FUN = ..parse_data,
    transformer = list(
      no_transformer,
      conv_transformer,
      rr_transformer,
      invar_transformer,
      invar_robust_transformer
    ),
    name = transformer_labels,
    MoreArgs = list("x" = x),
    SIMPLIFY = FALSE
  )
  names(data) <- transformer_labels

  # Parse transformed data
  transformed_data <- data.table::rbindlist(
    lapply(data, function(x) (x$transformed_data)),
    use.names = TRUE
  )
  transformed_data$transformation <- factor(
    transformed_data$transformation,
    levels = transformer_labels
  )

  # Parse residual data
  residual_data <- data.table::rbindlist(
    lapply(data, function(x) (x$residual_data)),
    use.names = TRUE
  )
  residual_data$transformation <- factor(
    residual_data$transformation,
    levels = transformer_labels
  )
  residual_data[, ":="("z_observed_truncated" = z_observed, "truncated" = FALSE)]
  residual_data[z_observed < plot_limits[1], ":="("z_observed_truncated" = plot_limits[1], "truncated" = TRUE)]
  residual_data[z_observed > plot_limits[2], ":="("z_observed_truncated" = plot_limits[2], "truncated" = TRUE)]

  # Parse central normality data.
  central_normality_data <- data.table::rbindlist(
    lapply(data, function(x) (x$central_normality_data)),
    use.names = TRUE
  )

  # Quantile-quantile plots.
  p_qq <- ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = .data$z_expected,
      y = .data$z_observed_truncated,
      colour = .data$transformation,
      size = .data$transformation,
      shape = .data$truncated
    ))
  p_qq <- p_qq + plot_theme
  p_qq <- p_qq + ggplot2::scale_colour_manual(
    name = "transformation",
    values = c("#111111", "#59a14f", "#8cd17d", "#e15759", "#ff9d9a"),
    breaks = transformer_labels,
    drop = FALSE
  )
  p_qq <- p_qq + ggplot2::scale_size_manual(
    name = "transformation",
    values = c(1.2, 1.2, 0.6, 1.2, 0.6),
    breaks = transformer_labels,
    drop = FALSE
  )
  p_qq <- p_qq + ggplot2::scale_shape_manual(
    values = c(16, 4),
    labels = c(FALSE, TRUE),
    guide = "none"
  )
  p_qq <- p_qq + ggplot2::coord_cartesian(
    xlim = plot_limits,
    ylim = plot_limits
  )
  p_qq <- p_qq + ggplot2::geom_abline(
    intercept = 0.0,
    slope = 1.0
  )
  p_qq <- p_qq + ggplot2::xlab("expected quantile")
  p_qq <- p_qq + ggplot2::ylab("observed quantile")

  if (!last_set) {
    p_qq <- p_qq + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      plot.margin = grid::unit(c(2, 2, 10, 2), "points")
    )
  }

  # Density plots.
  p_d <- ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = .data$x_standardised,
      y = ggplot2::after_stat(scaled),
      colour = .data$transformation
    ))
  p_d <- p_d + plot_theme
  p_d <- p_d + ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      margin = ggplot2::margin(0, 0, 2, 0)
    )
  )
  p_d <- p_d + ggplot2::theme(
    plot.margin = grid::unit(c(5, 2, 0, 2), "points")
  )

  p_d <- p_d + ggplot2::scale_colour_manual(
    name = "transformation",
    values = c("#111111", "#8cd17d", "#59a14f", "#ff9d9a", "#e15759"),
    breaks = transformer_labels,
    drop = FALSE,
    guide = "none"
  )
  p_d <- p_d + ggplot2::xlim(c(-3, 3))
  p_d <- p_d + ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )

  if (!first_set) {
    p_d <- p_d + ggplot2::theme(
      title = ggplot2::element_blank()
    )
  }

  # Without transformation -----------------------------------------------------
  p0_qq <- p_qq + ggplot2::geom_point(
    data = residual_data[transformation == "none"],
    show.legend = c(colour = TRUE, size = TRUE, shape = FALSE)
  )
  p0_d <- p_d + ggplot2::geom_density(
    data = transformed_data[transformation == "none"]
  )
  p0_d <- p0_d + ggplot2::labs(
    title = "original",
    tag = paste0(data_name)
  )

  # Standard transformation ----------------------------------------------------
  p1_qq <- p_qq + ggplot2::geom_point(
    data  = residual_data[transformation %in% c("conventional", "Raymaekers-Rousseeuw")],
    show.legend = c(colour = TRUE, size = TRUE, shape = FALSE)
  )
  p1_qq <- p1_qq + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )
  p1_d <- p_d + ggplot2::geom_density(
    data = transformed_data[transformation %in% c("conventional", "Raymaekers-Rousseeuw")]
  )
  p1_d <- p1_d + ggplot2::labs(title = "conventional")

  # Location- and scale invariant transformation -------------------------------
  p2_qq <- p_qq + ggplot2::geom_point(
    data  = residual_data[transformation %in% c("invariant", "robust invariant")],
    show.legend = c(colour = TRUE, size = TRUE, shape = FALSE)
  )
  p2_qq <- p2_qq + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )
  p2_d <- p_d + ggplot2::geom_density(
    data = transformed_data[transformation %in% c("invariant", "robust invariant")]
  )
  p2_d <- p2_d + ggplot2::labs(title = "invariant")

  return(list(
    "density" = list(p0_d, p1_d, p2_d),
    "qq" = list(p0_qq, p1_qq, p2_qq),
    "data" = data
  ))
}



.plot_marginal_effects <- function(plot_theme, manuscript_dir, subset = "numeric") {

  learner <- transformation_method <- normalisation_method <- task_difficulty <- NULL

  # Get ML experiment data.
  data <- .get_ml_experiment_data(manuscript_dir = manuscript_dir, subset = subset)

  # Drop single dataset classified as "unsolvable".
  data <- data[task_difficulty != "unsolvable"]
  data$task_difficulty <- droplevels(data$task_difficulty)

  # Build a random forest to create latent representations of the data.
  model_data <- data[, mget(c("learner", "transformation_method", "normalisation_method", "task_difficulty", "value_rank"))]
  model <- ranger::ranger(
    value_rank ~ learner + transformation_method + normalisation_method + task_difficulty,
    data = model_data,
    min.node.size = 2L,
    num.trees = 2000L,
    seed = 19L
  )

  # Generate all combinations of process parameters and task difficulty.
  abstract_data <- data.table::data.table(
    expand.grid(
      list(
        "learner" = levels(data$learner),
        "transformation_method" = levels(data$transformation_method),
        "normalisation_method" = levels(data$normalisation_method),
        "task_difficulty" = levels(data$task_difficulty)
      )
    )
  )
  abstract_data$task_difficulty <- factor(
    abstract_data$task_difficulty,
    levels=levels(abstract_data$task_difficulty),
    ordered=TRUE
  )

  # difficulty_levels <- c("all", "very_easy", "easy", "intermediate", "difficult", "very_difficult", "unsolvable")
  # difficulty_labels <- c("overall", "very easy", "easy", "intermediate", "difficult", "very difficult", "unsolvable")
  difficulty_levels <- c("all", "very_easy", "easy", "intermediate", "difficult", "very_difficult")
  difficulty_labels <- c("overall", "very easy", "easy", "intermediate", "difficult", "very difficult")

  # Effect of learner ----------------------------------------------------------
  result_list <- list()
  learners <- levels(data$learner)

  for (learner_a in learners) {
    for (learner_b in learners) {
      if (learner_a == learner_b) next

      x <- predict(model, data = abstract_data[learner == learner_a])$predictions -
        predict(model, data = abstract_data[learner == learner_b])$predictions

      result_list <- c(
        result_list,
        list(data.table::data.table(
          "value" = mean(x),
          "task_difficulty" = "all",
          "learner_1" = learner_a,
          "learner_2" = learner_b
        ))
      )

      for (difficulty in levels(abstract_data$task_difficulty)) {
        x <- predict(model, data = abstract_data[learner == learner_a & task_difficulty == difficulty])$predictions -
          predict(model, data = abstract_data[learner == learner_b & task_difficulty == difficulty])$predictions

        result_list <- c(
          result_list,
          list(data.table::data.table(
            "value" = mean(x),
            "task_difficulty" = difficulty,
            "learner_1" = learner_a,
            "learner_2" = learner_b
          ))
        )
      }
    }
  }

  data_1 <- data.table::rbindlist(result_list)
  data_1$task_difficulty <- factor(
    data_1$task_difficulty,
    levels = difficulty_levels,
    labels = difficulty_labels
  )
  data_1$learner_1 <- factor(
    data_1$learner_1,
    levels = learners
  )
  data_1$learner_2 <- factor(
    data_1$learner_2,
    levels = learners
  )


  # Effect of normalisation method ---------------------------------------------
  result_list <- list()
  x <- predict(model, data = abstract_data[normalisation_method == "robust_standardisation"])$predictions -
    predict(model, data = abstract_data[normalisation_method == "none"])$predictions
  result_list <- c(result_list, list(data.table::data.table("value" = mean(x), "task_difficulty" = "all")))

  for (difficulty in levels(abstract_data$task_difficulty)) {
    x <- predict(model, data = abstract_data[normalisation_method == "robust_standardisation" & task_difficulty == difficulty])$predictions -
      predict(model, data = abstract_data[normalisation_method == "none" & task_difficulty == difficulty])$predictions
    result_list <- c(result_list, list(data.table::data.table("value" = mean(x), "task_difficulty" = difficulty)))
  }
  data_2 <- data.table::rbindlist(result_list)
  data_2$task_difficulty <- factor(
    data_2$task_difficulty,
    levels = difficulty_levels,
    labels = difficulty_labels
  )


  # Effect of transformation method vs. none, by learner, task difficulty and
  # normalisation method.
  result_list <- list()
  for (lrnr in levels(abstract_data$learner)) {
    for (norm_method in levels(abstract_data$normalisation_method)){
      for (method in c("conventional", "invariant_robust", "invariant_robust_gof")) {
        x <- predict(model, data = abstract_data[learner == lrnr & transformation_method == method & normalisation_method == norm_method])$predictions -
          predict(model, data = abstract_data[learner == lrnr & transformation_method == "none" & normalisation_method == norm_method])$predictions
        x <- data.table::data.table(
          "value" = mean(x),
          "transformation_method" = paste0(method, "_none"),
          "normalisation_method" = norm_method,
          "learner" = lrnr,
          "task_difficulty" = "all"
        )
        result_list <- c(result_list, list(x))

        for (difficulty in levels(abstract_data$task_difficulty)) {
          x <- predict(model, data = abstract_data[learner == lrnr & task_difficulty == difficulty & transformation_method == method & normalisation_method == norm_method])$predictions -
            predict(model, data = abstract_data[learner == lrnr & task_difficulty == difficulty & transformation_method == "none" & normalisation_method == norm_method])$predictions
          x <- data.table::data.table(
            "value" = mean(x),
            "transformation_method" = paste0(method, "_none"),
            "normalisation_method" = norm_method,
            "learner" = lrnr,
            "task_difficulty" = difficulty
          )
          result_list <- c(result_list, list(x))
        }
      }

      for (method in c("invariant_robust", "invariant_robust_gof")) {
        x <- predict(model, data = abstract_data[learner == lrnr & transformation_method == method & normalisation_method == norm_method])$predictions -
          predict(model, data = abstract_data[learner == lrnr & transformation_method == "conventional" & normalisation_method == norm_method])$predictions
        x <- data.table::data.table(
          "value" = mean(x),
          "transformation_method" = paste0(method, "_conventional"),
          "normalisation_method" = norm_method,
          "learner" = lrnr,
          "task_difficulty" = "all"
        )
        result_list <- c(result_list, list(x))

        for (difficulty in levels(abstract_data$task_difficulty)) {
          x <- predict(model, data = abstract_data[learner == lrnr & task_difficulty == difficulty & transformation_method == method & normalisation_method == norm_method])$predictions -
            predict(model, data = abstract_data[learner == lrnr & task_difficulty == difficulty & transformation_method == "conventional" & normalisation_method == norm_method])$predictions
          x <- data.table::data.table(
            "value" = mean(x),
            "transformation_method" = paste0(method, "_conventional"),
            "normalisation_method" = norm_method,
            "learner" = lrnr,
            "task_difficulty" = difficulty
          )
          result_list <- c(result_list, list(x))
        }
      }
    }
  }
  data_3 <- data.table::rbindlist(result_list)
  data_3$learner <- factor(
    data_3$learner,
    levels = c("glm", "random_forest_ranger"),
    labels = c("GLM", "random forest")
  )
  transformation_method_labels = c(
    "conventional \nvs. none",
    "robust invariant \nvs. none",
    "robust invariant \nw. ECNT vs. none",
    "robust invariant \nvs. conventional",
    "robust invariant \nw. ECNT vs. conventional"
  )
  data_3$transformation_method <- factor(
    data_3$transformation_method,
    levels = c("conventional_none", "invariant_robust_none", "invariant_robust_gof_none", "invariant_robust_conventional", "invariant_robust_gof_conventional"),
    labels = transformation_method_labels
  )
  data_3$normalisation_method <- factor(
    data_3$normalisation_method,
    levels = c("none", "robust_standardisation"),
    labels = c("no normalisation", "z-standardisation")
  )
  data_3[, "learner_normalisation_method":= paste0(learner, "\n", normalisation_method)]
  data_3$learner_normalisation_method <- factor(
    data_3$learner_normalisation_method,
    levels = c("GLM\nno normalisation", "GLM\nz-standardisation", "random forest\nno normalisation", "random forest\nz-standardisation")
  )
  data_3$task_difficulty <- factor(
    data_3$task_difficulty,
    levels = difficulty_levels,
    labels = difficulty_labels
  )

  data_3_x_range = c(-0.08, 0.08)
  data_3_labels <- data.table::data.table(
    learner = factor("random forest", levels = levels(data_3$learner)),
    normalisation_method = factor("z-standardisation", levels = levels(data_3$normalisation_method)),
    transformation_method = factor(transformation_method_labels),
    task_difficulty  = factor("overall", levels = difficulty_labels),
    label_1 = c("none better", "none better", "none better", "conv. better", "conv. better"),
    label_2 = c("conv. better", "invar. better", "invar. better", "invar. better", "invar. better")
  )

  annotation_settings <- get_annotation_settings(plot_theme)

  p3 <- ggplot2::ggplot(
    data = data_3,
    mapping = ggplot2::aes(x = value, y = task_difficulty, fill = task_difficulty)
  )
  p3 <- p3 + plot_theme
  p3 <- p3 + ggplot2::theme(strip.clip = "off")
  p3 <- p3 + ggplot2::geom_vline(
    xintercept = 0.0,
    linetype = "longdash",
    colour = "grey40"
  )
  p3 <- p3 + ggplot2::geom_col(show.legend = FALSE)
  p3 <- p3 + ggplot2::facet_grid(learner_normalisation_method ~ transformation_method)
  p3 <- p3 + ggplot2::geom_text(
    x = data_3_x_range[1L],
    ggplot2::aes(label=label_1),
    hjust = "inward",
    vjust = "outward",
    data = data_3_labels,
    colour = annotation_settings$colour,
    family = annotation_settings$family,
    fontface = annotation_settings$face,
    size = annotation_settings$geom_text_size * 0.8
  )
  p3 <- p3 + ggplot2::geom_text(
    x = data_3_x_range[2L],
    ggplot2::aes(label=label_2),
    hjust = "inward",
    vjust = "outward",
    data = data_3_labels,
    colour = annotation_settings$colour,
    family = annotation_settings$family,
    fontface = annotation_settings$face,
    size = annotation_settings$geom_text_size * 0.8
  )
  p3 <- p3 + ggplot2::scale_fill_discrete(
    name = NULL,
    guide = "none",
    type = c(
      "unsolvable" = "#1c3319",
      "very difficult" = "#325d2d",
      "difficult" = "#498641",
      "intermediate" = "#60ae56",
      "easy" = "#87c280",
      "very easy" = "#add6a9",
      "overall" = "#e15759"
    )
  )
  p3 <- p3 + ggplot2::xlim(data_3_x_range)
  p3 <- p3 + ggplot2::ylab("task difficulty")
  p3 <- p3 + ggplot2::xlab("normalised rank difference")

  data_2_x_range = c(-0.005, 0.03)
  data_2_labels <- data.table::data.table(
    x = data_2_x_range,
    task_difficulty  = factor(c("overall", "overall"), levels = difficulty_labels),
    label = c("none better", "z-stand. better")
  )

  p2 <- ggplot2::ggplot(
    data = data_2,
    mapping = ggplot2::aes(x = value, y = task_difficulty, fill = task_difficulty)
  )
  p2 <- p2 + plot_theme
  p2 <- p2 + ggplot2::geom_vline(
    xintercept = 0.0,
    linetype = "longdash",
    colour = "grey40"
  )
  p2 <- p2 + ggplot2::geom_col(show.legend = FALSE)
  p2 <- p2 + ggplot2::geom_text(
    data = data_2_labels,
    mapping = ggplot2::aes(x = x, label = label),
    hjust = "inward",
    vjust = "outward",
    colour = annotation_settings$colour,
    family = annotation_settings$family,
    fontface = annotation_settings$face,
    size = annotation_settings$geom_text_size * 0.8
  )
  p2 <- p2 + ggplot2::scale_fill_discrete(
    name = NULL,
    guide = "none",
    type = c(
      "unsolvable" = "#482105",
      "very difficult" = "#813b08",
      "difficult" = "#bb550c",
      "intermediate" = "#f07014",
      "easy" = "#f4934e",
      "very easy" = "#f7b687",
      "overall" = "#e15759"
    )
  )
  p2 <- p2 + ggplot2::xlim(data_2_x_range)
  p2 <- p2 + ggplot2::ylab("task difficulty")
  p2 <- p2 + ggplot2::xlab("normalised rank difference")
  p2 <- p2 + ggplot2::ggtitle("z-standardisation vs. none")

  data_1_x_range = c(-0.05, 0.35)
  data_1_labels <- data.table::data.table(
    x = data_1_x_range,
    task_difficulty  = factor(c("overall", "overall"), levels = difficulty_labels),
    label = c("GLM better", "random forest better")
  )

  p1 <- ggplot2::ggplot(
    data = data_1,
    mapping = ggplot2::aes(x = value, y = task_difficulty, fill = task_difficulty)
  )
  p1 <- p1 + plot_theme
  p1 <- p1 + ggplot2::geom_vline(
    xintercept = 0.0,
    linetype = "longdash",
    colour = "grey40"
  )
  p1 <- p1 + ggplot2::geom_col(show.legend = FALSE)
  p1 <- p1 + ggplot2::geom_text(
    data = data_1_labels,
    mapping = ggplot2::aes(x = x, label = label),
    hjust = "inward",
    vjust = "outward",
    colour = annotation_settings$colour,
    family = annotation_settings$family,
    fontface = annotation_settings$face,
    size = annotation_settings$geom_text_size * 0.8
  )
  p1 <- p1 + ggplot2::scale_fill_discrete(
    name = NULL,
    guide = "none",
    type = c(
      "unsolvable" = "#192534",
      "very difficult" = "#2d435d",
      "difficult" = "#416186",
      "intermediate" = "#567fae",
      "easy" = "#809ec2",
      "very easy" = "#a9bed6",
      "overall" = "#e15759"
    )
  )
  p1 <- p1 + ggplot2::xlim(data_1_x_range)
  p1 <- p1 + ggplot2::ylab("task difficulty")
  p1 <- p1 + ggplot2::xlab("normalised rank difference")
  p1 <- p1 + ggplot2::ggtitle("random forest vs. generalised linear model")

  p <- (p1 + p2 + patchwork::plot_layout(ncol = 2, axes = "collect_y", axis_titles = "collect_y")) / p3 + patchwork::plot_layout(
    heights = c(0.3, 1.0)
  )

  return(p)
}


.plot_test_dependency_sample_size <- function() {
  n_rep <- 10L
  k <- 0.10

  # From 5 to 1000 in equal steps (log 10 scale)
  n <- 10^(seq(from = 0.75, to = 3.0, by = 0.125))

  set.seed(19L)

  data <- list()
  for (ii in seq_along(n)) {
    # Repeat 1000 times.
    for (jj in seq_len(n_rep)) {
      x <- power.transform::ragn(floor(n[ii]), location = 0, scale = 1 / sqrt(2), alpha = 0.5, beta = 2)

      # Compute upper and lower quartiles, and IQR.
      q_lower <- stats::quantile(x, probs = 0.25, names = FALSE)
      q_upper <- stats::quantile(x, probs = 0.75, names = FALSE)
      interquartile_range <- stats::IQR(x)

      # Set data where the outliers will be copied into.
      x_outlier <- x

      # Generate outlier values that are smaller than Q1 - 1.5 IQR or larger
      # than Q3 + 1.5 IQR.
      n_draw <- ceiling(k * length(x))
      x_random <- stats::runif(n_draw, min = -2.0, max = 2.0)


      outlier <- numeric(n_draw)
      if (any(x_random < 0)) {
        outlier[x_random < 0] <- q_lower - 1.5 * interquartile_range + x_random[x_random < 0] * interquartile_range
      }

      if (any(x_random >= 0)) {
        outlier[x_random >= 0] <- q_upper + 1.5 * interquartile_range + x_random[x_random >= 0] * interquartile_range
      }

      # Randomly insert outlier values.
      x_outlier[sample(seq_along(x), size = n_draw, replace = FALSE)] <- outlier

      data[[jj + (ii - 1L) * n_rep]] <- data.table::data.table(
        "n" = n[ii],
        "p_ecn" = power.transform::assess_transformation(
          x = x,
          transformer = power.transform::find_transformation_parameters(
            x = x,
            method = "none"
          ),
          verbose = FALSE
        ),
        "p_shap_wilk" = shapiro.test(x)$p.value,
        "p_ecn_outlier" = power.transform::assess_transformation(
          x = x_outlier,
          transformer = power.transform::find_transformation_parameters(
            x = x_outlier,
            method = "none"
          ),
          verbose = FALSE
        ),
        "p_shap_wilk_outlier" = shapiro.test(x_outlier)$p.value
      )
    }
  }

}
