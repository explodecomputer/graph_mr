library(igraph)


g <- make_ring(10)
distances(g)
shortest_paths(g, 5)
all_shortest_paths(g, 1, 6:8)
mean_distance(g)
## Weighted shortest paths
el <- matrix(nc=3, byrow=TRUE,
             c(1,2,0, 1,3,2, 1,4,1, 2,3,0, 2,5,5, 2,6,2, 3,2,1, 3,4,1,
               3,7,1, 4,3,0, 4,7,2, 5,6,2, 5,8,8, 6,3,2, 6,7,1, 6,9,1,
               6,10,3, 8,6,1, 8,9,1, 9,10,4) )
g2 <- add_edges(make_empty_graph(10), t(el[,1:2]), weight=el[,3])
distances(g2, mode="out")


adjm <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.9,0.1)), nc=10)
g1 <- graph_from_adjacency_matrix( adjm )

distances(g1)



mat <- matrix(rbind(
	c(0,1,0,1,0,0),
	c(0,0,1,0,0,0),
	c(0,0,0,0,0,1),
	c(0,0,0,0,1,0),
	c(0,0,0,0,0,1),
	c(0,0,0,0,0,0)), 6,6
)

g1 <- graph_from_adjacency_matrix(mat, mode="directed")
distances(g1)
g1

all_shortest_paths(g1, 1, 6)



mat <- matrix(rbind(
	c(0,1,0,0,0,0),
	c(0,0,1,0,0,0),
	c(0,0,0,1,0,0),
	c(0,0,0,0,1,0),
	c(0,0,0,0,0,1),
	c(0,0,0,0,0,0)), 6,6
)

g1 <- graph_from_adjacency_matrix(mat, mode="directed")
distances(g1)
g1

all_shortest_paths(g1, from=1, to=6)

mat <- matrix(runif(36), 6, 6)
g1 <- graph_from_adjacency_matrix(mat < 0.15, mode="directed")
res <- all_shortest_paths(g1, 1, 6)


