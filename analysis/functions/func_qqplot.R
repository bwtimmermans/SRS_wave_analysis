# qqplot function derived from extRemes package.
# Modified to allow changing font size of legend and line size of regression fit.

func_qqplot <- function (x, y, pch = 19, xlab = "x Quantiles", ylab = "y Quantiles", 
    regress = TRUE, make.plot = TRUE, cex.leg = 3.0, lwd.reg = 3.0, ...) 
{
    args <- list(...)
    out <- list()
    out$call <- match.call()
    out$names <- list(x = as.character(substitute(x)), y = as.character(substitute(y)))
    x <- sort(na.omit(x))
    y <- sort(na.omit(y))
    qy <- extRemes::quantilefun(y)
    m <- length(x)
    n <- length(y)
    N <- m + n
    M <- m * n/N
    K <- 1.36
    p <- (1:m - 1)/(m - 1)
    yq <- qy(p)
#print(paste("qy:",qy(p - K/sqrt(M))))
    yl <- qy(p - K/sqrt(M))
    yu <- qy(p + K/sqrt(M))
    if (make.plot) {
        if (is.null(args$xlim) && is.null(args$ylim)) 
            plot(x, yq, pch = pch, xlim = range(x), ylim = range(yq, 
                yl, yu, na.rm = TRUE), xlab = xlab, ylab = ylab, 
                ...)
        else if (is.null(args$xlim)) 
            plot(x, yq, pch = pch, xlim = range(x), xlab = xlab, 
                ylab = ylab, ...)
        else if (is.null(args$ylim)) 
            plot(x, yq, pch = pch, ylim = range(yq, yl, yu, na.rm = TRUE), 
                xlab = xlab, ylab = ylab, ...)
        else plot(x, yq, pch = pch, xlab = xlab, ylab = ylab, 
            ...)
        lines(x, yl, lty = 2, col = "gray", lwd = lwd.reg)
        lines(x, yu, lty = 2, col = "gray", lwd = lwd.reg)
        abline(0, 1, lty = 2, col = "darkorange", lwd = lwd.reg)
    }
    if (regress) {
        fit <- lm(y ~ x, data = data.frame(x = x, y = yq))
        if (make.plot) {
            lines(x, predict(fit), col = "grey", lty = 1, lwd = lwd.reg)
            #legend("topleft", legend = c("1-1 line", "regression line", 
            #    "95% confidence bands"), col = c("darkorange", 
            #    "grey", "gray"), lty = c(2, 1, 2), bty = "n", cex = cex.leg)
            legend("topleft", legend = c("1-1 line", "regression line","95% confidence bands"),
		   col = c("darkorange","grey","grey"), lty = c(2, 1, 2), bty = "n", cex = cex.leg)
        }
        out$regression <- fit
    }
    else if (make.plot) 
        legend("bottomright", legend = c("1-1 line", "95% confidence bands"), 
            col = c("darkorange", "gray"), lty = c(2, 2), bty = "n", cex = cex.leg)
    out$qdata <- data.frame(x = x, y = yq, lower = yl, upper = yu)
    class(out) <- "qqplot"
    invisible(out)
}

