 # MixMPLN
"MixMPLN" is a package written in R, which has two features. First feature is generating a proper synthetic sample-taxa count matrix. Second feature is receiving a sample-taxa count matrix and extracting k(number of components) different interaction networks between taxa.  

## 1-Prerequisites
Following packages must be installed and loaded in R environment before using "MixMPLN":
```
library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(cvTools)
library(huge)
library(edgeR)
library(factoextra)
library(igraph)
library(ROCR)
```
## 2-Installation
After installing the previous prerequisities, please run the following lines to install the "MixMPLN"
```
library("devtools")
install_github("sahatava/MixMPLN")
library("MixMPLN")
```
## 3-Generating Synthetic data
"MixMPLN" is able to generate the synthetic matrix which follows the real data features. If you are interested to apply the "MixMPLN" on your own data you could skip this part and jump to the next part (4-Extracting K interaction networks). 
Our generated synthetic data consists of K componenets. For each componenet, a random positive definite matrix with specific level of sparsity as a real precision matrix is generated. N generated samples for each componenet follows the MPLN distribution based on the mentioned precision matrix and a random mean vector. Finally, samples from different componets are combined together to generate multi component original count data.
The generated original count data follows by multinomial distribution to mimic the sampling process.   
##### Usage
```
out_generate = generate(K , N , d , sp, type)
```
##### Arguments
```
K-------> number of components
N-------> number of samples in each component
d-------> number of taxa for the synthetic sample-taxa count matrix
sp------> a value between 0 and 1 (sparsity level)
type----> data type: "orig" (no sampling applied) , "samp" (after sampling) , "TMM" (after normalization)
```

##### Values
```
out_generate$M -----------------> synthetic sample-taxa count matrix(rows represent samples and columns represent taxa)
out_generate$real_precision[[i]]-----> i_th real precision matrix 
```
## 4-Extracting K interaction networks
MixMPLN is a fuction which accept a sample-taxa matrix as an input to extract K different underlied interarction networks. Please save your data as matrix which its rows represents samples and its columns represents the taxa names. If you are interested in visualising the final graphs, you need to save the taxa names as colnames of this input matrix. 
K( number of components ) is another input for MixMPLN which should be specified by the user.   
For initialize the optimization using the Kmeans clustering please specified init="Kmeans" and rep=1. Otherwise, init="Random" and rep can be any integer number. As higher value for rep results in higher running time, we suggest to set it as 3 to trade off between the time and accuracy.  
##### Usage
```
out_MixMPLN = MixMPLN(M , K , penalty , init , rep )
```
##### Arguments
```
M--------> name of the csv file which includes the sample-taxa matrix
K-----------> number of components 
penalty------> method to select tuning parameter: "no_sparsity" , "CV" (cross validation) , "StARS" , "fixed" , "iterative" 
init---------> initialization can be "Kmeans" or "Random"
rep----------> number of repeats when initialization is "Random"
out---------> the type of output for each componenet: precision matrix, partial correlation, lambda, mixing coefficiant, clustering membership 
threshold---> a value between 0 and 1. This threshold is used to generate adjacency matrix from partial correlation
```

##### Values
```
out_MixMPLN$precision[[i]]-----> i_th precision matrix(dXd)
out_MixMPLN$partial[[i]]-------> i_th partial correlation matrix(dXd)
out_MixMPLN$lambda[[i]]--------> i_th lambda matrix(NXd)
out_MixMPLN$pi-----------------> mixing coefficiant
out_MixMPLN$cluster------------> clustering membership
```

## 5-Auxiliary functions

### silhouette method to find the optimal K
```
fviz_nbclust(M , MixMPLN ,penalty ,k.max, method = c("silhouette"))
```
##### Arguments
```
M-----------> sample-taxa count matrix
K.max-------> k=1 to k=k.max will be checked
penalty------> method to select tuning parameter: "no_sparsity" , "CV" (cross validation) , "StARS" , "fixed" , "iterative" 
```
##### Values
```
Running the mentioned command will plot a diagram which shows optimum K 
```
 
### Graphical display
```
visualize(partial , threshold )
```
##### Arguments
```
partial--------> partial correlation matrix corresponding to one component
threshold------> Threshold to get the adjacency matrix of the graph
```
##### Values
```
The graphical plot for the taxa interaction network. 
```

### Converting the partial corelation to the adjacency matrix
```
out_adj = partial2adj(partial, threshold)
```
##### Arguments
```
partial-------------> partial correlation matrix as a result of running the MixMPLN
threshold-----------> a value between 0 and 1
```
##### Values
```
out_adj-------------> 1 and -1 indicate edge existance, 0 indicates no edge 
```

### Comparing the estimated precision matrix and the real one for synthetic data
```
out = compare(out_generate$real_precision , out_MixMPLN$precision)
```
##### Arguments
```
ut_generate$real_precision -------> real partial correlation
out_MixMPLN$precision -------> estimated partial correlation
```
##### Values
```
out$frob------> forbenious norm of the difference
out$fr--------> relative difference
out$rms-------> rms difference
out$sn--------> sensitivity
out$sp--------> specificity
out$ROC-------> area under the ROC 
```


## 6-Example
### Example on synthetic data
```
out_generate = generate(K=2 , N=30 , d=30 , sp=0.8, type="orig")
out_MixMPLN = MixMPLN( out_generate$M , K=2 , penalty="CV" , init = "Random" , rep = 2)
out = compare(out_generate$real_precision  , out_MixMPLN$precision)
``` 
### Example on real data
```
M=as.matrix(read.csv("real_data.csv",sep=",",header=TRUE,row.names=1,check.names=FALSE))
fviz_nbclust(M , MixMPLN ,"CV",k.max, method = c("silhouette"))
out_MixMPLN = MixMPLN( M , K=2 , penalty="CV", init = "Random" , rep = 2)
visualize(out_MixMPLN$partial[[1]] , threshold=0.3 )
```

## Publication
The "MixMPLN" method is explained in following paper:
Sahar Tavakoli and Shibu Yooseph. Learning a mixture of microbial networks using
minorization-maximization. To appear in Intelligent Systems for Molecular Biology, 2019
Please note that in our implementation we have used the "huge" package in all the variations. 
