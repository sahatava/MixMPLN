
library(igraph)
 


######################
#' A function to    
#'    
#' @param input
#' @param ts
#' @return res
#' @export
partial2adj <- function(input, ts){
res = input
for(i in 1:nrow(input)){
  for(j in 1:ncol(input)){
    if(input[i,j]> ts){res[i,j]=1}
    if(input[i,j]< -ts){res[i,j]=-1}
    if(input[i,j]>= -ts && input[i,j]<= ts ){res[i,j]=0}
  }
}
return(res)
}
######################
#' A function to    
#'    
#' @param input
#' @param ts
#' @return res
#' @export
pos <- function(input, ts){
res = input
for(i in 1:nrow(input)){
  for(j in 1:ncol(input)){
    if(input[i,j]> ts){
      res[i,j]=1
    }else{
      res[i,j]=0
    }
    
  }
}
return(res)
}
################
#' A function to    
#'    
#' @param input
#' @param ts
#' @return res
#' @export
neg <- function(input, ts){
res = input
for(i in 1:nrow(input)){
  for(j in 1:ncol(input)){
    if(input[i,j] < -ts){
      res[i,j]=1
    }else{
      res[i,j]=0
    }
    
  }
}
return(res)
}
################
#' A function to    
#'    
#' @param partial
#' @param treshold
#' @return res
#' @export

visualize <- function(partial , threshold  ){

a1 = partial2adj(partial , threshold)
a1pos = pos(partial , threshold)
a1neg = neg(partial , threshold)

colnames(a1) = colnames(partial) 
 
colnames(a1pos) = colnames(partial)
colnames(a1neg) = colnames(partial) 
 
 
g1pos = graph_from_adjacency_matrix(a1pos, mode = "undirected", weighted=T, diag = FALSE, add.colnames = NULL, add.rownames = NA)
g1neg = graph_from_adjacency_matrix(a1neg, mode = "undirected", weighted=T, diag = FALSE, add.colnames = NULL, add.rownames = NA)

E(g1pos)$color = 'green'
E(g1neg)$color = 'magenta'


tree_union = graph.union(g1pos, g1neg, byname=T)
TU_col = edge_attr(tree_union, "color_1")
E2 = which(is.na(edge_attr(tree_union, "color_1")))
TU_col[E2] = edge_attr(tree_union, "color_2")[E2]
g1new = set_edge_attr(tree_union, "color", value=TU_col)
 

g = g1new
### Example
## Generate some fake data
#n <- 75
#g <- erdos.renyi.game(n, 0.5)
#V(g)$name = paste("long_name", 1:n, sep="_")

## laid out as a circle to begin with
la <- layout.circle(g)

par(mar=c(10,10,10,10))
plot(g, layout=la, vertex.size=1, vertex.label="")

## Apply labels manually
#Specify x and y coordinates of labels, adjust outward as desired
x = la[,1]*1.55
y = la[,2]*1.55

#create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
angle = ifelse(atan(-(la[,1]/la[,2]))*(180/pi) < 0,  90 + atan(-(la[,1]/la[,2]))*(180/pi), 270 + atan(-la[,1]/la[,2])*(180/pi))

#Apply the text labels with a loop with angle as srt
for (i in 1:length(x)) {
	text(x=x[i], y=y[i], labels=V(g)$name[i], adj=NULL, pos=NULL, cex=.8, col="black", srt=angle[i], xpd=T)
        #print(i)
        #print(V(g)$name[i])
}

return()

}


