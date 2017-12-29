library("ape")
tr <- read.tree(file="/Users/yuanenming/Desktop/phylogenetic-tree-construction/test.out")
plot(tr, main="phylogenetic Tree", tip.color = hsv(runif(15, 0.65, 0.95), 1, 1, 0.7), edge.color = hsv(runif(10, 0.65, 0.75), 1, 1, 0.7), edge.width = runif(20, 0.5, 3), use.edge.length = TRUE, col = "gray80")
axisPhylo()
