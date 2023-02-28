.check_ggtheme <- function(ggtheme) {

  # Check if the provided theme is a suitable theme.
  if (inherits(ggtheme, "theme")) {

    # Check if the theme is complete.
    if (!attr(ggtheme, "complete")) {
      stop(paste0(
        "The plotting theme is not complete. The most likely cause is lack of a valid template, such as ",
        "theme_familiar or ggplot2::theme_light. Note that ggplot2::theme is designed to tweak existing ",
        "themes when creating a plot."))
    }

  } else if (is.null(ggtheme)) {
    ggtheme <- ggplot2::theme_light(base_size = 9)

  } else {
    stop(paste0("The provided ggplot2 theme is not a theme. Found: ", class(ggtheme)))
  }

  return(ggtheme)
}
