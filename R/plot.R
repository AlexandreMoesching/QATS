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

#' Displays the current fit along with the data
#'
#' @param xx.0 True hidden path
#' @param xx.1 1st estimated path
#' @param xx.2 2nd 2stimated path
#' @param par Model parameters
#' @param SS Partition
#' @param k_new Newly added change points
#'
#' @return Plot
#' @keywords internal
display.0 <- function(xx.0, xx.1 = NULL, xx.2 = NULL,
                      par, SS = NULL, k_new = NULL) {
  # Load some variables
  mseq <- par$mseq
  n <- par$n

  # Prepare plot window
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 0.2, 0.5))
  plot(1:n, xx.0,
    pch = 16, yaxt = "n", cex = 1,
    ylim = range(mseq),
    xlab = expression(italic(k)),
    ylab = expression(italic(x[k]))
  )
  graphics::axis(side = 2, labels = mseq, at = mseq)
  if (!is.null(xx.2)) graphics::lines(1:n, xx.2, col = "red", lwd = 4)
  if (!is.null(xx.1)) graphics::lines(1:n, xx.1, col = "blue", lwd = 2)
  # Partition
  if (!is.null(SS)) graphics::abline(v = SS[, 1], lty = 2)
  # New change points
  if (!is.null(k_new)) graphics::abline(v = k_new, lwd = 2, col = "green")

  # Data
  plot(1:n, par$yy,
    pch = 16, cex = 0.5,
    xlab = expression(italic(k)),
    ylab = expression(italic(y[k]))
  )

  # There are no other plots...
}

#' Displays the current fit and the gain functions, along with the data
#'
#' @param xx.0 True hidden path
#' @param xx.1 1st estimated path
#' @param xx.2 2nd 2stimated path
#' @param par Model parameters
#' @param SS Partition
#' @param k_new Newly added change points
#' @param CP True change points
#' @param l Left-most index of the current window
#' @param res1 Vector of 1-dimensional gain
#' @param t_star1 True maximum
#' @param k_star1 Maximum found
#' @param res2 Matrix of 2-dimensional gain
#' @param t_star2 True maxium
#' @param k_star2 Maximum found
#'
#' @return Plot
#' @keywords internal
display.1 <- function(xx.0, xx.1 = NULL, xx.2 = NULL,
                      par, SS = NULL, k_new = NULL,
                      CP = NULL, l = 1,
                      res1 = NULL, t_star1 = NULL, k_star1 = NULL,
                      res2 = NULL, t_star2 = NULL, k_star2 = NULL) {
  # Load some variables
  mseq <- par$mseq
  n <- par$n

  # Prepare plot window
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 0.2, 0.5))
  plot(1:n, xx.0,
    pch = 16, yaxt = "n", cex = 1,
    ylim = range(mseq),
    xlab = expression(italic(k)),
    ylab = expression(italic(x[k]))
  )
  graphics::axis(side = 2, labels = mseq, at = mseq)
  if (!is.null(xx.2)) graphics::lines(1:n, xx.2, col = "red", lwd = 4)
  if (!is.null(xx.1)) graphics::lines(1:n, xx.1, col = "blue", lwd = 2)
  # Partition
  if (!is.null(SS)) graphics::abline(v = SS[, 1], lty = 2)
  # New change points
  if (!is.null(k_new)) graphics::abline(v = k_new, lwd = 2, col = "green")

  # Data
  plot(1:n, par$yy,
    pch = 16, cex = 0.5,
    xlab = expression(italic(k)),
    ylab = expression(italic(y[k]))
  )

  # 1D function to maximize
  if (!is.null(res1) && length(unique(c(res1[!is.na(res1)]))) >= 1) {
    plot((l - 1) + (1:length(res1)), res1,
      type = "l",
      xlab = expression(italic(k[1])), ylab = "Score"
    )
    graphics::abline(v = CP)
    # True maximum
    if (!is.null(t_star1)) {
      graphics::abline(v = (l - 1) + t_star1, lwd = 2, col = "black")
    }
    # Maximum found
    if (!is.null(k_star1)) {
      # cat("Maximum found (red)\n",
      #   "Press [enter] to continue\n",
      #   sep = ""
      # )
      graphics::abline(v = k_star1, col = "red")
    }
  }

  # 2D function to maximize
  if (!is.null(res2) && length(unique(c(res2[!is.na(res2)]))) > 1) {
    plot3D::image2D(res2,
      col = grDevices::hcl.colors(100, "Oslo"), # hcl.pals()
      # col = grDevices::topo.colors(100),
      x = (l - 1) + (1:nrow(res2)),
      y = (l - 1) + (1:ncol(res2)),
      xlab = expression(italic(k[1])),
      ylab = expression(italic(k[2]))
    )
    graphics::abline(h = CP)
    graphics::abline(v = CP)
    # True maximum
    if (!is.null(t_star2)) {
      graphics::points((l - 1) + t_star2, col = "black", pch = 16, cex = 1)
    }
    # Maximum found
    if (!is.null(k_star2)) {
      # cat("Maximum found (red)\n",
      #   "Press [enter] to continue\n",
      #   sep = ""
      # )
      graphics::points(k_star2, col = "red", pch = 16, cex = 1)
    }
  }
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
    graphics::plot(l:r, res,
      type = "l",
      xlab = expression(italic(k[1])),
      ylab = "Gain"
    )
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
    plot3D::image2D(res,
      # col = grDevices::topo.colors(100),
      col = grDevices::hcl.colors(100, "Oslo"), # hcl.pals()
      # contour = TRUE,
      x = l:r,
      y = l:r,
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
#' @return Plot
#' @export
display.result <- function(xx.0, xx.1 = NULL, xx.2 = NULL,
                           par, yy = NULL) {
  kk <- which(diff(xx.0) != 0)
  ii <- c(rbind(c(1, kk + 1), c(kk, par$n)))
  if (!is.null(yy)) {
    myylim <- range(yy)
    myylab <- expression(italic(list(x[k], y[k])))
  } else {
    myylim <- range(par$mseq)
    myylab <- expression(italic(x[k]))
  }
  plot(ii, xx.0[ii],
    type = "n",
    yaxt = "n", ylim = myylim,
    xlab = expression(italic(k)),
    ylab = myylab
  )
  if (!is.null(yy)) {
    graphics::points(yy, col = "grey", cex = 0.2)
  }
  graphics::lines(ii, xx.0[ii], type = "l", col = "black", lwd = 4)
  graphics::axis(side = 2, labels = par$mseq, at = par$mseq)
  if (!is.null(xx.2)) {
    kk <- which(diff(xx.2) != 0)
    ii.2 <- c(rbind(c(1, kk + 1), c(kk, par$n)))
    graphics::lines(ii.2, xx.2[ii.2], col = "red", lwd = 2)
  }
  if (!is.null(xx.1)) {
    kk <- which(diff(xx.1) != 0)
    ii.1 <- c(rbind(c(1, kk + 1), c(kk, par$n)))
    graphics::lines(ii.1, xx.1[ii.1], col = "blue", lwd = 1)
  }
}
