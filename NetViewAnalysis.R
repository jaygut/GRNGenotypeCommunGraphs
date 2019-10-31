library(netview)
library(DT)
library(tidyverse)
library(igraph)
library(cowplot)

baseDir <- "/Users/jaysongutierrez/Desktop/MyPortfolioGitRepos/GRNGenotypeCommunGraphs/"
setwd(baseDir)

#---------------------------------------------------------------------------
# Custom functions ...

communityGraph <- function(outputGraph,plotTitle,communLayout, vertexSize){
  #'Plot graph highlighting communities.
  #'By default, the community layout is set to mds
  #'communLayout options: layout.mds, layout.fruchterman.reingold, layout.kamada.kawai
  coords <- layout.auto(outputGraph)
  
  plot(outputGraph, 
       vertex.size=vertexSize, 
       vertex.label=NA,
       mark.shape = 0.75,
       mark.groups= cluster_infomap(outputGraph),
       mark.border=NA,
       layout=communLayout(outputGraph),
       main = plotTitle,
       rescale=TRUE,
       xlim=c(-1,1), 
       ylim=c(-1,1))
}

run_Netview <- function(distMatrix,mkNNG,nodeColorVect) {
  #'This function runs the NetView pipeline on a genetic distance matrix 
  #'and returns the computed community graph and a df with some individual properties
  #'Inputs: 
  #'distMatrix: genetic distance 
  #'mkNNG: number of k neighbors
  #'nodeColorVect: vector of colors to be assigned to nodes 
  # Set num individuals for formatting metadata
  
  numInds = dim(distMatrix)[1];
  
  # Set options for graph
  popOptions <- netviewOptions(selectionTitle="NA", 
                               nodeID="ID", nodeGroup="Group", nodeColour="Colour", 
                               communityAlgorithms=c("Infomap"))
  
  # Metadata df for the algorithm: colnames = "ID", "Group", "Colour"
  popMetaData <- data.frame("ID" = as.character(seq(1,numInds)),
                            "Group" = rep(paste0("DerivedPop",derivedPop,"AncBackground",ancBackground),numInds),
                            "Colour" = nodeColorVect #rep("white",numInds)
  )
  
  # Compute graph based on genetic distance 
  graphs <- netview::netview(distMatrix, popMetaData, k=mkNNG, cluster = TRUE, options=popOptions)
  
  getGraph <- graphs[[1]] 
  
  #Calculate local centrality measure, degree centrality, which counts the number of 
  #links held by each node and points at individuals who can quickly connect with the 
  #wider network
  degreeCent = centr_degree(getGraph)$res
  
  #Calculate average nearest neighbor degree
  avgNND <- knn(getGraph)$knn
  
  #Calculate the clustering coefficient per vertex 
  clusteringCoeff <- transitivity(getGraph, type = "local")
  
  #Get community membership per node
  commun.membership <- get.graph.attribute(getGraph)$infomap$membership
  
  # Cast graph properties into a df
  df <- data.frame(IndIndex=1:numInds, CommunIndex=commun.membership)
  df <- df[with(df, order(IndIndex)),]
  row.names(df) <- NULL
  df$DegreeCentrality <- degreeCent
  df$avgNearestNeighDegree <- avgNND
  df$clusteringCoefficient <- clusteringCoeff
  
  output <- list(outputGraph = getGraph, graphPropDF = df)
  
}


getNodeColorVect4DuplInds <- function(duplIndFeats, ancBackground, derivedPop, numInds){
  #'Use this function to construct a list of colors for the graph based on info for
  #'the individuals duplicated for the gene in our in silico GRN evolution experiments
  
  # Let's construct nodeColorVect so that one can visualize in the population graph the individuals whose
  # genotype was initially duplicated for the Activator gene. In this way one has to considered 27 individuals
  # to be colored differently than the rest given that we ran 27 simulation replicates with randomly picked individuals
  # duplicated.
  indDuplIDs <- duplIndFeats[duplIndFeats$anc_gc==as.numeric(ancBackground) & duplIndFeats$derived_pop==as.numeric(derivedPop),
                             c("dup_carrying_indid")]
  
  # Also color code the nodes representing initially duplicated genotypes with red if the duplication event in that
  # simulation run failed to undergo fixation, or green if succeded to fix!
  #First, let's check in which of the cases the duplicate remained fixed >= 50 generations
  longestDuplRemainedFixed <- duplIndFeats[duplIndFeats$anc_gc==as.numeric(ancBackground) & duplIndFeats$derived_pop==as.numeric(derivedPop),
                                           c("longest_time_window_dupl_nearly_fixed")]
  
  indDuplColorCode <- c()
  for (n in 1:length(longestDuplRemainedFixed)) {
    if (longestDuplRemainedFixed[n]>=50) {
      indDuplColorCode[n] <- "green"
    } else {
      indDuplColorCode[n] <- "red"
    }
  }
  
  nodeColorVect <- c()
  for (i in 1:numInds) {
    if (i %in% indDuplIDs) {
      #color code individuals initially duplicated for the Act gene as red (unfixed) or green (fixed)
      nodeColorVect[i] <- indDuplColorCode[match(i,indDuplIDs)]
    } else{
      #non-duplicated genotypes are white color-coded
      nodeColorVect[i] <- "white"
    }
    
  }
  
  return(nodeColorVect)
}

#---------------------------------------------------------------------------

# Setting input file name and loading it in
derivedPop <- "1";
ancBackground <- "1";
targetFileName <- stringr::str_interp("GeneticDistanceMatrix_DerivedPop${derivedPop}_FromAnc${ancBackground}.txt")
folderPath <- "./popgen/"
fileNamePath <- paste0(folderPath,targetFileName)
outputFN <- stringr::str_interp("PopGeneticStructGraphs_DerivedPop${derivedPop}_FromAnc${ancBackground}.pdf")
plotTitle <- stringr::str_interp("Derived Pop. ${derivedPop}")
# Load in file containing the optimal knn for each instance. This was computed in Python using the package kneed, which
# compute the point in axis x at which max curvature is observed in the curve of knn vs num of communities
optimal_knn <- read_csv('./optimal_knn.csv')
opt.knn_vector <- optimal_knn[[as.numeric(ancBackground)]]
mkNNG <- opt.knn_vector[as.numeric(derivedPop)] #num of knn
vertexSize <- 3

# Load in file containing some pre-computed features for the initially duplicated genotypes in each simulation replicate 
duplIndFeatsFN <- paste0(baseDir,"./InitDuplGeneFixationFeats.csv")
duplIndFeats <- read.table(duplIndFeatsFN,header = TRUE, sep = ",")

#------------------------------------------------------------------------

# Run example for one particular population
#nodeColorVect <- rep("white",numInds)
distMatrix <- as.matrix(read.table(fileNamePath,sep = ","))
numInds <- dim(distMatrix)[1]

nodeColorVect <- getNodeColorVect4DuplInds(duplIndFeats, ancBackground, derivedPop, numInds)

graphProps <- run_Netview(distMatrix,mkNNG,nodeColorVect)

#save individual community graph plot as pdf
plotTitle = " "
graphPlotFN <- stringr::str_interp("./PopGeneticStructGraph4DerivedPop${derivedPop}FromAnc${ancBackground}.pdf")
pdf(paste0(baseDir,graphPlotFN)) 
communityGraph(graphProps$outputGraph, 
               plotTitle, 
               layout_with_kk,#layout.fruchterman.reingold,
               vertexSize)
dev.off()
