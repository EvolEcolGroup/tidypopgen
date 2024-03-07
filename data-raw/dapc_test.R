x <- glSim(50,4e3-50, 50, ploidy=2)
x
plot(x)

## perform DAPC
dapc1 <- dapc(x, n.pca=10, n.da=1)
dapc1

## plot results
scatter(dapc1, scree.da=FALSE)

## SNP contributions
loadingplot(dapc1$var.contr)
loadingplot(tail(dapc1$var.contr, 100), main="Loading plot - last 100 SNPs")

test_gt <- as_gen_tibble(x)
x <- test_gt
dapc2 <- dapc(x, n.pca=10, n.da=1)
dapc2
scatter(dapc2, scree.da=FALSE)
## SNP contributions
loadingplot(dapc2$var.contr)
loadingplot(tail(dapc2$var.contr, 100), main="Loading plot - last 100 SNPs")
