function (coo, xlim, ylim, border = "#333333", col = NA, lwd = 1, 
    lty = 1, points = FALSE, first.point = TRUE, cex.first.point = 0.5, 
    centroid = TRUE, xy.axis = TRUE, pch = 1, cex = 0.5, main = NA, 
    poly = TRUE, plot.new = TRUE, plot = TRUE, zoom = 1, ...) 
{
    coo <- coo_check(coo)
    if (plot.new) {
        op <- par(mar = c(3, 3, 2, 1))
        on.exit(par(op))
        if (!missing(zoom)) {
            half.side <- max(apply(coo, 2, function(x) diff(range(x))))/2
            add.xy <- half.side * zoom
            orig.xy <- coo_centpos(coo)
            xlim <- c(orig.xy[1] - add.xy, orig.xy[1] + add.xy)
            ylim <- c(orig.xy[2] - add.xy, orig.xy[2] + add.xy)
        }
        if (!missing(xlim) | !missing(ylim)) {
            if (missing(xlim)) {
                xlim <- ylim
            }
            if (missing(ylim)) {
                ylim <- xlim
            }
            plot(coo, type = "n", asp = 1, las = 1, cex.axis = 2/3, 
                ann = FALSE, frame = FALSE, xlim = xlim, ylim = ylim)
        }
        else {
            plot(coo, type = "n", asp = 1, las = 1, cex.axis = 2/3, 
                ann = FALSE, frame = FALSE)
        }
        if (xy.axis) {
            abline(h = 0, v = 0, col = "grey80", lty = 2)
        }
    }
    if (plot) {
        if (!missing(poly)) {
            if ((!poly) & missing(points)) 
                points <- TRUE
        }
        if (poly) {
            polygon(coo, col = col, border = NA)
            lines(coo, col = border, lwd = lwd, lty = lty)
        }
        if (missing(points)) {
            if (nrow(coo) <= 60) 
                points <- TRUE
        }
        if (points) {
            points(coo, pch = pch, cex = cex, col = border)
        }
        if (first.point) {
            angle <- atan2(coo[2, 2] - coo[1, 2], coo[2, 1] - 
                coo[1, 1]) * (180/pi) - 90
            text(coo[1, 1], coo[1, 2], labels = "^", cex = cex.first.point, 
                srt = angle)
        }
        if (centroid) {
            cent <- coo_centpos(coo)
            points(cent[1], cent[2], pch = 3, col = border, cex = cex)
        }
        if (!missing(main)) 
            title(main = main)
    }
}
