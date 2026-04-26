#' Gain function for path with 1 segment
#'
#' @param D Vector of indices
#' @param s Length of D
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Scalar representing the gain function on D
#' @keywords internal
H1_dbl <- function(D, s, x0, par) {
  res <- H1_dbl_cpp(
    D[1] - 1, D[s] - 1, x0 - 1, par$m,
    par$logPi, par$qq, par$GG
  )
  return(res)
}

#' Gain function for path with 2 segments
#'
#' @param D Vector of indices
#' @param s Length of D
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Vector representing the gain function on D
#' @keywords internal
H2_vec <- function(D, s, x0, par) {
  res <- H2_vec_cpp(
    D[1] - 1, D[s] - 1, x0 - 1, par$m,
    par$logPi, par$qq, par$GG
  )
  t_star <- which(res == max(res, na.rm = TRUE))
  return(list(res = res, t_star = t_star))
}

#' Gain function for path with 3 segments
#'
#' @param D Vector of indices
#' @param s Length of D
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Matrix representing the gain function on D
#' @keywords internal
H3_mat <- function(D, s, x0, par) {
  res <- H3_mat_cpp(
    D[1] - 1, D[s] - 1, x0 - 1, par$m,
    par$logPi, par$qq, par$GG
  )
  t_star <- which(res == max(res, na.rm = TRUE), arr.ind = TRUE)
  return(list(res = res, t_star = t_star))
}

# ── Helper: build step-function segments from a state vector ──────────────────

.step_segments <- function(xx, n) {
  kk <- which(diff(xx) != 0)
  l  <- c(1L, kk + 1L)
  r  <- c(kk, n)
  data.frame(l = l, r = r, x = xx[l])
}

#' Displays the current fit along with the data
#'
#' @param xx.0 True hidden path
#' @param xx.1 1st estimated path
#' @param xx.2 2nd estimated path
#' @param par Model parameters
#' @param SS Partition
#' @param k_new Newly added change points
#'
#' @return Plot
#' @keywords internal
display.0 <- function(xx.0, xx.1 = NULL, xx.2 = NULL,
                      par, SS = NULL, k_new = NULL) {
  n    <- par$n
  mseq <- par$mseq

  # Panel 1: state sequences
  seg0 <- .step_segments(xx.0, n)
  p1   <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg0,
      ggplot2::aes(x = l, xend = r, y = x, yend = x),
      linewidth = 1
    ) +
    ggplot2::scale_y_continuous(breaks = mseq) +
    ggplot2::labs(x = expression(italic(k)), y = expression(italic(x[k])))

  if (!is.null(xx.2)) {
    seg2 <- .step_segments(xx.2, n)
    p1 <- p1 + ggplot2::geom_segment(
      data = seg2,
      ggplot2::aes(x = l, xend = r, y = x, yend = x),
      colour = "#e34a33", linewidth = 0.8
    )
  }
  if (!is.null(xx.1)) {
    seg1 <- .step_segments(xx.1, n)
    p1 <- p1 + ggplot2::geom_segment(
      data = seg1,
      ggplot2::aes(x = l, xend = r, y = x, yend = x),
      colour = "#3182bd", linewidth = 0.5
    )
  }
  if (!is.null(SS))
    p1 <- p1 + ggplot2::geom_vline(
      xintercept = SS[, 1], linetype = "dashed", colour = "grey50"
    )
  if (!is.null(k_new))
    p1 <- p1 + ggplot2::geom_vline(
      xintercept = k_new, linewidth = 0.8, colour = "#31a354"
    )

  # Panel 2: observations
  df_y <- data.frame(k = seq_len(n), y = par$yy)
  p2   <- ggplot2::ggplot(df_y, ggplot2::aes(k, y)) +
    ggplot2::geom_point(size = 0.3, colour = "grey50") +
    ggplot2::labs(x = expression(italic(k)), y = expression(italic(y[k])))

  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 0.2, 0.5))
  print(p1)
  print(p2)
}

#' Displays the current fit and the gain functions, along with the data
#'
#' @param xx.0 True hidden path
#' @param xx.1 1st estimated path
#' @param xx.2 2nd estimated path
#' @param par Model parameters
#' @param SS Partition
#' @param k_new Newly added change points
#' @param CP True change points
#' @param l Left-most index of the current window
#' @param res1 Vector of 1-dimensional gain
#' @param t_star1 True maximum
#' @param k_star1 Maximum found
#' @param res2 Matrix of 2-dimensional gain
#' @param t_star2 True maximum
#' @param k_star2 Maximum found
#'
#' @return Plot
#' @keywords internal
display.1 <- function(xx.0, xx.1 = NULL, xx.2 = NULL,
                      par, SS = NULL, k_new = NULL,
                      CP = NULL, l = 1,
                      res1 = NULL, t_star1 = NULL, k_star1 = NULL,
                      res2 = NULL, t_star2 = NULL, k_star2 = NULL) {
  n    <- par$n
  mseq <- par$mseq

  # Panel 1: state sequences
  seg0 <- .step_segments(xx.0, n)
  p1   <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg0,
      ggplot2::aes(x = l, xend = r, y = x, yend = x),
      linewidth = 1
    ) +
    ggplot2::scale_y_continuous(breaks = mseq) +
    ggplot2::labs(x = expression(italic(k)), y = expression(italic(x[k])))

  if (!is.null(xx.2)) {
    seg2 <- .step_segments(xx.2, n)
    p1 <- p1 + ggplot2::geom_segment(
      data = seg2,
      ggplot2::aes(x = l, xend = r, y = x, yend = x),
      colour = "#e34a33", linewidth = 0.8
    )
  }
  if (!is.null(xx.1)) {
    seg1 <- .step_segments(xx.1, n)
    p1 <- p1 + ggplot2::geom_segment(
      data = seg1,
      ggplot2::aes(x = l, xend = r, y = x, yend = x),
      colour = "#3182bd", linewidth = 0.5
    )
  }
  if (!is.null(SS))
    p1 <- p1 + ggplot2::geom_vline(
      xintercept = SS[, 1], linetype = "dashed", colour = "grey50"
    )
  if (!is.null(k_new))
    p1 <- p1 + ggplot2::geom_vline(
      xintercept = k_new, linewidth = 0.8, colour = "#31a354"
    )

  # Panel 2: observations
  df_y <- data.frame(k = seq_len(n), y = par$yy)
  p2   <- ggplot2::ggplot(df_y, ggplot2::aes(k, y)) +
    ggplot2::geom_point(size = 0.3, colour = "grey50") +
    ggplot2::labs(x = expression(italic(k)), y = expression(italic(y[k])))

  # Panel 3: 1D gain
  p3 <- NULL
  if (!is.null(res1) && length(unique(c(res1[!is.na(res1)]))) >= 1) {
    df1 <- data.frame(k1 = (l - 1) + seq_along(res1), score = res1)
    p3  <- ggplot2::ggplot(df1, ggplot2::aes(k1, score)) +
      ggplot2::geom_line() +
      ggplot2::labs(x = expression(italic(k[1])), y = "Score")

    if (!is.null(CP))
      p3 <- p3 + ggplot2::geom_vline(xintercept = CP, colour = "grey50")
    if (!is.null(t_star1))
      p3 <- p3 + ggplot2::geom_vline(
        xintercept = (l - 1) + t_star1, linewidth = 0.8
      )
    if (!is.null(k_star1))
      p3 <- p3 + ggplot2::geom_vline(
        xintercept = k_star1, colour = "#e34a33"
      )
  }

  # Panel 4: 2D gain
  p4 <- NULL
  if (!is.null(res2) && length(unique(c(res2[!is.na(res2)]))) > 1) {
    if (!requireNamespace("plot3D", quietly = TRUE)) {
      stop("Package 'plot3D' is required for 2D gain plots. ",
           "Install it with install.packages('plot3D').")
    }
    # plot3D::image2D does not have a ggplot2 equivalent that is as
    # compact, so we keep it as a base-graphics panel for this one case.
    plot3D::image2D(res2,
      col = grDevices::hcl.colors(100, "Oslo"),
      x   = (l - 1) + seq_len(nrow(res2)),
      y   = (l - 1) + seq_len(ncol(res2)),
      xlab = expression(italic(k[1])),
      ylab = expression(italic(k[2]))
    )
    if (!is.null(CP)) {
      graphics::abline(h = CP)
      graphics::abline(v = CP)
    }
    if (!is.null(t_star2))
      graphics::points((l - 1) + t_star2, pch = 16, cex = 1)
    if (!is.null(k_star2))
      graphics::points(k_star2, col = "#e34a33", pch = 16, cex = 1)
  }

  # Arrange panels
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 0.2, 0.5))
  print(p1)
  print(p2)
  if (!is.null(p3)) print(p3)
  # p4 is printed by plot3D directly
}

#' Plot of the gain function for path with 2 segments
#'
#' @param l Left-most index
#' @param r Right-most index
#' @param x0 Previous state
#' @param par Model parameters
#' @param display Whether or not to display the map
#'
#' @return Plot
#' @export
display.vec <- function(l, r, x0, par, display = TRUE) {
  res <- H2_vec_cpp(
    l - 1, r - 1, x0 - 1,
    par$m, par$logPi, par$qq, par$GG
  )
  if (display) {
    df <- data.frame(k1 = l:r, gain = res)
    p  <- ggplot2::ggplot(df, ggplot2::aes(k1, gain)) +
      ggplot2::geom_line() +
      ggplot2::labs(x = expression(italic(k[1])), y = "Gain")
    print(p)
  }
  return(res)
}

#' Plot of the gain function for path with 3 segments
#'
#' @param l Left-most index
#' @param r Right-most index
#' @param x0 Previous state
#' @param par Model parameters
#' @param include.dim1 If true, add the single change point score
#' @param include.dim0 If true, add the no change point score
#' @param display Whether or not to display the map
#'
#' @return Plot
#' @export
display.mat <- function(l, r, x0, par,
                        include.dim1 = FALSE,
                        include.dim0 = FALSE,
                        display = TRUE) {
  res <- H3_mat_cpp(
    l - 1, r - 1, x0 - 1,
    par$m, par$logPi, par$qq, par$GG
  )
  if (include.dim1) {
    tmp <- H2_vec_cpp(
      l - 1, r - 1, x0 - 1,
      par$m, par$logPi, par$qq, par$GG
    )
    res[1, ] <- tmp
  }
  if (include.dim0) {
    diag(res) <- H1_dbl_cpp(
      l - 1, r - 1, x0 - 1,
      par$m, par$logPi, par$qq, par$GG
    )
  }
  if (display) {
    if (!requireNamespace("plot3D", quietly = TRUE)) {
      stop("Package 'plot3D' is required for 2D gain plots. ",
           "Install it with install.packages('plot3D').")
    }
    plot3D::image2D(res,
      col = grDevices::hcl.colors(100, "Oslo"),
      x   = l:r,
      y   = l:r,
      xlab = expression(italic(k[1])),
      ylab = expression(italic(k[2]))
    )
  }
  return(res)
}

#' Display final results
#'
#' @param xx.0 True hidden path
#' @param xx.1 1st estimated path
#' @param xx.2 2nd estimated path
#' @param par Model parameters
#' @param yy The observation sequence
#'
#' @return A ggplot object (also printed)
#' @export
display.result <- function(xx.0, xx.1 = NULL, xx.2 = NULL,
                           par, yy = NULL) {
  n    <- par$n
  mseq <- par$mseq
  seg0 <- .step_segments(xx.0, n)

  if (!is.null(yy)) {
    myylim <- range(yy)
    myylab <- expression(italic(x[k] * " / " * y[k]))
  } else {
    myylim <- range(mseq)
    myylab <- expression(italic(x[k]))
  }

  p <- ggplot2::ggplot() +
    ggplot2::scale_y_continuous(breaks = mseq, limits = myylim) +
    ggplot2::labs(x = expression(italic(k)), y = myylab)

  if (!is.null(yy)) {
    df_y <- data.frame(k = seq_len(n), y = yy)
    p <- p + ggplot2::geom_point(
      data = df_y, ggplot2::aes(k, y),
      colour = "grey70", size = 0.2
    )
  }

  # True path (black, thick)
  p <- p + ggplot2::geom_segment(
    data = seg0,
    ggplot2::aes(x = l, xend = r, y = x, yend = x),
    linewidth = 1.2
  )

  if (!is.null(xx.2)) {
    seg2 <- .step_segments(xx.2, n)
    p <- p + ggplot2::geom_segment(
      data = seg2,
      ggplot2::aes(x = l, xend = r, y = x, yend = x),
      colour = "#e34a33", linewidth = 0.7
    )
  }
  if (!is.null(xx.1)) {
    seg1 <- .step_segments(xx.1, n)
    p <- p + ggplot2::geom_segment(
      data = seg1,
      ggplot2::aes(x = l, xend = r, y = x, yend = x),
      colour = "#3182bd", linewidth = 0.4
    )
  }

  print(p)
  invisible(p)
}
