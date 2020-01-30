# define live cell gate based on 'openCyto' 
# plug-in from http://opencyto.org/plugins.html
# x = expression matrix, q = quantile, bs = line intercept & slope
#' @importFrom flowCore polygonGate
#' @importFrom grDevices chull
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats qnorm
#' @importFrom openCyto register_plugins
.live_gate <- function(x, q = 0.99, bs = c(1, 0.5)) {
    # specifying gating function
    .gating_fun <- function(fr, pp_res, channels = NA, id = "", ...) {
        # subset channels of interest
        x <- exprs(fr[, channels]) 
        # scale data for comparison w/ 'qnorm()'
        x0 <- scale(x) 
        # set boundary level as q-th quantile of standard normal
        z <- qnorm(q) 
        # find p(x) for that level
        pd <- dmvnorm(c(z, z))[1] 
        px <- dmvnorm(x0)
        # find points above boundary level 
        keep1 <- px > pd 
        # find points below line y = a + b * x  
        keep2 <- (bs[1] + bs[2] * x0[, 1]) > x0[, 2] 
        # intersection of points below line & above threshold level
        pts <- x[keep1 & keep2, ] 
        # get boundary points (convex hull) 
        pts <- pts[chull(pts), ]
        # return gate
        polygonGate(.gate = pts, filterId = id)
    }
    # register gate
    suppressMessages(
        foo <- register_plugins(
            fun = .gating_fun, 
            methodName = "liveGate", 
            dep = "mvtnorm", 
            "gating"))
}
.get_elip <- function(mu, cov, q, n = 250) {
    es <- eigen(cov)
    e <- es$vec %*% diag(sqrt(es$val))
    r <- sqrt(qchisq(q, df = 2))
    theta <- seq(0, 2*pi, l = n)
    v <- cbind(r * cos(theta), r * sin(theta))
    xy <- data.frame(t(mu - (e %*% t(v))))
    colnames(xy) <- c("x", "y"); xy
}
#' @importFrom dplyr bind_rows do group_by group_modify rename_at
#' @importFrom purrr map
.get_gate <- function(gs, type, group_by, ...) {
    args <- list(...)
    switch(type, 
        rect = {
            df <- data.frame(
                do.call(rbind, map(gs, "min")),
                do.call(rbind, map(gs, "max")))
            colnames(df) <- c("xmin", "ymin", "xmax", "ymax")
            if (!is.null(group_by)) 
                df[[group_by]] <- names(gs)
            return(df)
        },
        elip = {
            es <- lapply(gs, function(u) .get_elip(u@mean, u@cov, args$q))
            bind_rows(map(es, data.frame), .id = group_by)
        },
        live = {
            xy <- map(gs, "boundaries")
            df <- map(xy, data.frame)
            df <- bind_rows(df, .id = group_by)
            chs <- setdiff(colnames(df), group_by)
            df <- rename_at(df, chs, ~c("x", "y"))
            # add 1st point for each group
            # to avoid disconnected path
            if (is.null(group_by)) {
                df <- rbind(df, head(df, n = 1))
            } else {
                df <- group_by(df, !!sym(group_by)) %>% 
                    group_modify(~head(.x, 1L)) %>% 
                    do(rbind(df, .data))
            }
            return(df)
        }
    )

}
