x <- glSim(100,500,500)
x
plot(x)
grp <- find.clusters(x, stat = "BIC") # choose =FALSE
plot(grp$Kstat, type = "o", xlab = "number of clusters (K)",
     ylab = "BIC",
     main = "find.clusters on a genlight object\n(two groups)")
