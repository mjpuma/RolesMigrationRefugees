# Set working directory.
setwd("D:/Documents/Work/University of Florida/Postdoc 2021-2022/Project Documents/Roles of Mobility Networks")

# turning off scientific notation
options(scipen = 999)

# function to normalize from 0 to 1 for plotting arguments
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# set custom colors
custom_colors <- c("#9e0142", 
                   "#d53e4f",
                   "#f46d43",
                   "#fdae61",
                   "#fee08b",
                   "#e6f598",
                   "#abdda4",
                   "#66c2a5",
                   "#3288bd",
                   "#5e4fa2")

# Make sure to have the following packages installed before using library() function:
library(readxl)
library(igraph) # Remember that igraph conflicts with statnet & sna.
library(sna)
library(blockmodeling)
library(FactoMineR)
library(factoextra)
library(NbClust)
# library(RColorBrewer)
# library(MASS)
# library(vegan)
library(ggplot2)
library(tidyverse)
library(ggraph) # need to manually install tidygraph. Select "no" for tidygraph installation when asked to install from source requiring compilation
# library(stats)
library(gridExtra)
# library(fpc) # select "no" for kernlab installation when asked to install from source requiring compilation
# library(scales)
library(clValid)
library(xUCINET)
library(mclogit)
library(car)
library(ggpubr)
# library(ggplotify)

### country reference from bilat_mig
States <- read.csv("data/country_list.csv", header = T)

# Migrant Flow Estimates 1990-2020
edges.mig <- read.csv("data/bilat_mig.csv", header = T)[,c(1:3,9)]
edges.mig$state.origin.name <- States$name[match(edges.mig$orig, States$iso3c)]
edges.mig$state.destination.name <- States$name[match(edges.mig$dest, States$iso3c)]
edges.mig <- edges.mig[,c(1,5,6,4)]
edges.mig$year0 <- edges.mig$year0+5
colnames(edges.mig)[1] <- "year5"
edges.mig$from <- paste0(edges.mig$state.origin.name, ".", edges.mig$year5)
edges.mig$to <- paste0(edges.mig$state.destination.name, ".", edges.mig$year5)
edges.mig <- edges.mig[,c(5,6,4)]

##### Migrant networks - strict cutoff - unweighted #####
a <- na.omit(edges.mig[edges.mig[,3]>=1,])
colnames(a)[3] <- "weight"
a1 <- graph_from_data_frame(a)
a2 <- as.matrix(as_adjacency_matrix(a1, type = "both", names = T))
a3 <- graph_from_adjacency_matrix(a2)
a4 <- igraph::delete_vertices(a3, V(a3)[igraph::degree(a3)==0])
a5 <- as.matrix(as_adjacency_matrix(a4, type = "both", names = T))
a6 <- graph_from_adjacency_matrix(a5)
a7 <- as.matrix(as_adj(a6))
# a8 <- redist(a7, mode = "digraph")
# a8 <- xREGE(a7)
a8 <- REGE.nm.for(a7, E = 1)
a9 <- a8$E
a9 <- 1-a9
# rownames(a8) <- colnames(a8) <- rownames(a7)

# block <- blockmodel(a7, ec = equiv.clust(dat = a7, equiv.dist = a8, mode = "digraph", diag = FALSE, cluster.method = "average"), h=0.99)
# plot(block)
link <- "average"
set.seed(123)
plot(hcut(x = a9, hc_method = link))

### hierarchical clustering assessment
k_medoids_plots <- list()
# elbow method
k_medoids_plots[[1]] <- fviz_nbclust(a9,
                                     FUNcluster = hcut,
                                     hc_method = link,
                                     diss = a9,
                                     nstart = 1,
                                     k.max = 25,
                                     method = "wss") +
  labs(title = "UNDESA", subtitle = "Elbow method")

# silhouette method
k_medoids_plots[[2]] <- fviz_nbclust(a9,
                                     FUNcluster = hcut,
                                     hc_method = link,
                                     diss = a9,
                                     nstart = 1,
                                     k.max = 25,
                                     method = "silhouette") +
  labs(subtitle = "Silhouette method")

# gap statistic
# set.seed(123)
# k_medoids_plots[[3]] <- fviz_nbclust(a9,
#                                      FUNcluster = hcut,
#                                      diss = as.dist(a9),
#                                      method = "gap_stat",
#                                      nstart = 1,
#                                      k.max = 50,
#                                      nboot = 500) +
#   labs(subtitle = "Gap statistic method")

do.call(grid.arrange, c(k_medoids_plots, ncol = 1))

for (i in 2:5){
  # Hierarchical Clustering
  p <- fviz_cluster(hcut(x = a9, hc_method = link, k = i), data = a9) +
    labs(title = "UNDESA", subtitle = paste0(link," link: k = ",i, sep="")) +
    theme(legend.position="none")
  # K-medoids
  # p <- fviz_cluster(pam(x = a9, diss = T, k = i), data = a9) + 
  #   labs(title = "UNDESA", subtitle = paste0("k-medoids: k = ",i, sep="")) +
  #   theme(legend.position="none")
  plot(p)
}


a <- hcut(x = a9, method = link)
# a$height <- round(a$height, 10)
b <- cutree(a, k = 3)
b
c <- data.frame("id"=names(b),"role"=b)


write.csv(c,"results/UNDESA_full_roles_k3.csv")
clusters <- read.csv("results/UNDESA_full_roles_k3.csv")[,2:3]

### Network Metrics
mig_network_stats <- data.frame()
mig_mean_network_stats <- data.frame()

a <- na.omit(edges.mig[edges.mig[,3]>=1,])
colnames(a)[3] <- "weight"
a1 <- graph_from_data_frame(a)
a2 <- as.matrix(as_adjacency_matrix(a1, type = "both", names = T, attr = "weight"))
a3 <- graph_from_adjacency_matrix(a2, weighted = T)
a4 <- igraph::delete_vertices(a3, V(a3)[igraph::degree(a3)==0])
a5 <- as.matrix(as_adjacency_matrix(a4, type = "both", names = T, attr = "weight"))
a6 <- graph_from_adjacency_matrix(a5, weighted = T)
in.degree <- data.frame(in.degree = igraph::degree(a6, mode = "in"))
out.degree <- igraph::degree(a6, mode = "out")
in.strength <- igraph::strength(a6, mode = "in")
out.strength <- igraph::strength(a6, mode = "out")
av.in.strength <- igraph::strength(a6, mode = "in")/igraph::degree(a6, mode = "in")
av.out.strength <- igraph::strength(a6, mode = "out")/igraph::degree(a6, mode = "out")
betweenness <- igraph::betweenness(a6, directed = TRUE, weights = E(a6)$weight, nobigint = FALSE)
region <- rownames(in.degree)
mig_network_stats <- cbind(region,in.degree,out.degree,in.strength,out.strength,av.in.strength,av.out.strength,betweenness)
mig_network_stats <- merge(mig_network_stats, clusters, by.x = "region", by.y = "id")
b <- list()
for (j in 1:length(unique(na.omit(mig_network_stats[,9])))){
  b[[j]] <- c(j,
              nrow(mig_network_stats[mig_network_stats[,9]==j,]),
              mean(mig_network_stats[mig_network_stats[,9]==j,]$in.degree),
              mean(mig_network_stats[mig_network_stats[,9]==j,]$out.degree),
              mean(mig_network_stats[mig_network_stats[,9]==j,]$in.strength),
              mean(mig_network_stats[mig_network_stats[,9]==j,]$out.strength),
              mean(mig_network_stats[mig_network_stats[,9]==j,]$av.in.strength),
              mean(mig_network_stats[mig_network_stats[,9]==j,]$av.out.strength),
              mean(mig_network_stats[mig_network_stats[,9]==j,]$betweenness))
}
mig_mean_network_stats <- do.call("rbind", b)
colnames(mig_mean_network_stats) <- c("role","n","in.degree","out.degree","in.strength","out.strength","av.in.strength","av.out.strength","betweenness")
rm(a,av.in.strength,av.out.strength,a1,a2,a3,a4,a5,a6,b,betweenness,in.degree,in.strength,out.degree,out.strength,region)

write.csv(mig_mean_network_stats,"results/UNDESA_full_roles_k3_metrics.csv")
