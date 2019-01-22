setwd("C:/Users/neo204/OneDrive - University of Exeter/NEVO/Rushton Model/Paper/Tidy Code/")

rm(list = ls())

outputs = read.table("plot_model_runs.txt")

par(mfrow = c(1,1), family = "sans", mar = c(4,4,1,1), ps = 14)

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 3000), ylim = c(0, 30))
axis(side = 1, at = c(0, 1000, 2000, 3000), labels = TRUE)
mtext("time [day]", side = 1, line = 2.5)
axis(side = 2, at = c(0, 5, 10, 15, 20, 25, 30), labels = TRUE)
mtext(expression("river flow [" ~ m^3/s ~ "]"), side = 2, line = 2.5)
lines(1:2922, outputs[1,], col = "red")
lines(1:2922, outputs[2,], col = "orange")
lines(1:2922, outputs[3,], col = "yellow")
lines(1:2922, outputs[4,], col = "green")
lines(1:2922, outputs[5,], col = "blue")
legend("topleft", legend = c(expression((list(x[1], x[2], x[3], x[4]))),
                             expression((list(1, 0, 0, 0))),
                             expression((list(1/4, 1/4, 1/4, 1/4))),
                             expression((list(0, 1/3, 1/3, 1/3))),
                             expression((list(1/2, 0, 0, 1/2))),
                             expression((list(0, 0, 0, 1)))),
       col = c("black", "red", "orange", "yellow", "green", "blue"),
       lty = c(0, rep(1, 5)), lwd = rep(2, 6), bty = "n")
box(lwd = 1)