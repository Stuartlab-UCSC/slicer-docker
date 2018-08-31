# Run slicer on the example data.
library(SLICER)
library("lle")
data(traj)
genes = select_genes(traj)
cell_names <- colnames(traj)
if (is.null(cell_names)) {
  cell_names <- 1:dim(traj)[1]
}
k = select_k(traj[,genes], kmin=5)
traj_lle = lle(traj[,genes], m=2, k)$Y

traj_graph = conn_knn_graph(traj_lle,5)
ends = find_extreme_cells(traj_graph, traj_lle)
start = 1
cell_ordering = cell_order(traj_graph, start)
branch_assignment = assign_branches(traj_graph,start)
save(branch_assignment, cell_ordering, cell_names, file="/home/rstudio/shared/example_obj.rda")
####################################################3
# Make the common json format
# graph structure is made by assuming all branches are
# connected by pseudotime.
nodes <- c()
g <- igraph::make_empty_graph(directed=F)
stem_name <- "stem"
nodes <- c(nodes, stem_name)

# Makes a star graph with node names
g <- igraph::add_vertices(g,1, name=stem_name)
for (branch in unique(branches)){
  node_name <- paste(c("end_", branch), collapse = "")
  nodes <- c(nodes, stem_name)
  g <- igraph::add_vertices(g, 1, name=node_name)
  g <- igraph::add_edges(g, c(stem_name, node_name))
}

edges <- igraph::as_edgelist(g)
edgeIds <- c()
nodeIds1 <- c()
nodeIds2 <- c()
# CM stands for cell mapping.
cellIdCM <- c()
edgeIdCM <- c()
pseudotimeCM <- c()
for (row in 1:(dim(edges)[1])){
  node1 <- edges[row,1]
  nodeIds1 <- c(nodeIds1, node1)
  node2 <- edges[row,2]
  nodeIds2 <- c(nodeIds2, node1)
  edge_id <- paste(c(node1,"_", node2), collapse="")
  edgeIds <- c(edgeIds, edge_id)
  branch_n <- as.numeric(tail(strsplit(node2, "-branch-")[[1]],1))
  n_branch_cells <- sum(branches==branch_n)
  edgeIdCM <- c(edgeIdCM, rep(edge_id, n_branch_cells))
  cellIdCM <- c(cellIdCM, cell_names[branches==branch_n])
  pseudotimeCM <- c(pseudotimeCM, cells_ordered[branches==branch_n])
}
output <- list(
  nodes= list(nodeId=nodes),
  egdes= list(
    edgeId=edgeIds,
    nodeId1= edges[,1],
    nodeId2= edges[,2]
  ),
  cellMapping=list(
    cellId= cellIdCM,
    edgeId=edgeIdCM,
    psuedotime=pseudotimeCM
  )
)


############################################################
# Look the pseudotime is linear how really does that work
# with the branches?
plot(g)
library(car)
scatterplot(x=traj_lle[,1], y=traj_lle[,2],labels=branches)
colnames(traj_lle) <- c('x', "y")
dim(traj_lle)
length(branches)
?data.frame
df = data.frame(traj_lle[,1], traj_lle[,2], branches)
colnames(df) <- c('x', "y", "branches")
library(ggplot2)
ggplot2.scatterplot(
  data=df,
)

ggplot(df, aes(x=x, y=y, color=branches)) + geom_point(shape=1)
ggplot(df, aes(x=x, y=y, color=cells_ordered)) + geom_point(shape=1)
order_check <- cells_ordered
# Check for the median and color all cells greater than it as 1
greater_than_2_median <- rep(0,length(order_check))
greater_than_2_median[order_check > (median(order_check[branches == 2])+30)] <- 1
ggplot(df, aes(x=x, y=y, color=greater_than_2_median)) + geom_point(shape=1)

# Could progmatically test for that with correlation of
# distance matrix of the x-y points with the lowest psuedo time in a branch,
# have to decide a cutoff, OR if there's 3 branches can take the lowest and assume
# that  it's the middle branch.
b2<-traj_lle[branches==2]
minb2ind <- which(cells_ordered==min(b2), arr.ind=T)
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
distToMin2 <- function(x) {
  return(euc.dist(traj_lle[minb2ind], x))
}
distance_from_first_smallest <- sapply(b2,distToMin2)
max(distance_from_first_smallest)
which()



##############################################
source("/home/rstudio/shared/slicer_convert.r")
a <- to_cell_x_branch(branch_assignment, cell_ordering, 1:500)
b <- to_common_list(branch_assignment, cell_ordering, cell_names)
