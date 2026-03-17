

plot_expression_density_kde_ks <- function(df,
                                           x_col = "x",
                                           y_col = "y",
                                           expr_col = "expr",
                                           gridsize = c(300, 300),
                                           xmin = NULL, xmax = NULL,
                                           ymin = NULL, ymax = NULL,
                                           H = NULL,
                                           adjust = 1,
                                           transform = c("none", "log1p"),
                                           guide_title = "Weighted KDE") {
  transform <- match.arg(transform)
  stopifnot(all(c(x_col, y_col, expr_col) %in% names(df)))
  
  x <- df[[x_col]]
  y <- df[[y_col]]
  w <- pmax(df[[expr_col]], 0)
  
  ok <- is.finite(x) & is.finite(y) & is.finite(w)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  
  if (is.null(xmin)) xmin <- min(x)
  if (is.null(xmax)) xmax <- max(x)
  if (is.null(ymin)) ymin <- min(y)
  if (is.null(ymax)) ymax <- max(y)
  mins <- c(xmin, ymin)
  maxs <- c(xmax, ymax)
  
  X <- cbind(x, y)
  
  if (is.null(H)) {
    H <- Hpi(X, binned = TRUE) * adjust
  } else {
    H <- H * adjust
  }
  
  kd <- kde(x = X,
            w = w,
            H = H,
            gridsize = gridsize,
            xmin = mins,
            xmax = maxs,
            binned = TRUE)
  
  gx <- kd$eval.points[[1]]
  gy <- kd$eval.points[[2]]
  z  <- kd$estimate
  
  den_df <- expand.grid(x = gx, y = gy)
  den_df$value <- as.vector(z)
  
  if (transform == "log1p") {
    den_df$value <- log1p(den_df$value)
  }
  
  ggplot(den_df, aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_viridis_c(
      name = guide_title,
      option = "magma"
    ) +
    theme_void() +
    theme(legend.position = "right")
}



df <- data.frame(
  x = dat_meta$y, 
  y = max(dat_meta$x)-dat_meta$x,
  expr = c(dat_meta$NECTIN3)
)

plot_expression_density_kde_ks(df,
                               gridsize = c(400, 400),
                               transform = "none",
                               guide_title = "exp")+labs(subtitle = "NECTIN3")
