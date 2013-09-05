#################################################################################
# Runner script for flowCL which tests functionality -- this would be the 
# script to convert into an R package eventually
# Author: Radina Droumeva (radina.droumeva@gmail.com)
# Date last modified: July 19, 2013
# Edited by: Justin Meskas (jmeskas@bccrc.ca)
# Date last modified: Sept 05, 2013
#
# The idea is to start with a phenotype, such as "CD4+" or "CD4+CD8+", and get 
# back the relevant portions of the ontology tree.
# There is a preliminary scoring system in place which should be assessed and 
# changed appropriately to aid in selecting the best cell label match.

# 
# Clean memory
rm(list = ls())

# Variable used for computer programming efficiency calculating
tmpf = tempfile()
Rprof(tmpf)

initialTime <- Sys.time()

VisualizationSkip       <- TRUE     # Skips the making of the Visualization (tree diagram)
penalty.skip            <- TRUE     # Skips the penalty calculation (slowest part of the code, and not that important)
ShowOntologyMarkerNames <- FALSE    # Shows either the common marker name or the ontology marker name in the tree diagram

# Load all phenotypes
 listPhenotypes <- read.csv("CL_Markers/CL_Markers_Singles.csv", header=F)            # Individual Markers
#  listPhenotypes <- read.csv("CL_Markers/CL_Markers_HIPC_Multiples.csv", header=F)     # Multiple Markers
# listPhenotypes <- read.csv("CL_Markers/CL_Markers_HIPC_Everything.csv", header=F)    # Individual Markers and Multiple Markers
# listPhenotypes <- read.csv("CL_Markers/CL_Markers_AllMultiples", header=F)           # Non HIPC Combinations of markers
# listPhenotypes <- read.csv("CL_Markers/CL_Markers_Aug_27_2013_Test.csv", header=F)   # List of Markers to create fig for presentation Sept/03/2013  
listPhenotypes <- as.matrix(listPhenotypes)

# Start and end of the iterations through the "CL_Markers_().csv" file
IterStart = 1
# IterEnd   = 6
IterEnd   = length(listPhenotypes)

fname <- "CL_Markers/markers_NickName_OntologyName.csv"
if(!file.exists(fname)){
    temp.table <- matrix()
    temp.table[1] <-paste('CD4+');    temp.table[2] <-paste('CD4 molecule')
#     temp.table <- t(temp.table)
    write.table(t(temp.table), fname, sep = ",", col.names = FALSE, row.names = FALSE)  
}

# Preallocate lists
listPhenotypeUpdate  <- listPhenotypes
listPhenotypeOriginal<- listPhenotypes
listPhenotypeID      <- listPhenotypes
listCellID           <- listPhenotypes
listPhenotypeSuccess <- listPhenotypes
listMarkerLabels     <- listPhenotypes
listCellLabels       <- listPhenotypes

for(q in IterStart:IterEnd){    
  listPhenotypeSuccess[[q]] <-   "No"
  listPhenotypeID[[q]]      <-   "Nothing"
  listCellID[[q]]           <-   "Nothing"
  listMarkerLabels[[q]]     <-   "Nothing"
  listCellLabels[[q]]       <-   "Nothing"
}
                          
# Define directories for storing results:
save.dir             <- 'results/'
save.dirResults      <- 'results/results/'
save.dirParents      <- 'results/parents/'
save.dirParentsQuery <- 'results/parents_query/'

dir.create(save.dir,             showWarnings=FALSE, recursive=TRUE)
dir.create(save.dirResults,      showWarnings=FALSE, recursive=TRUE)
dir.create(save.dirParents,      showWarnings=FALSE, recursive=TRUE)
dir.create(save.dirParentsQuery, showWarnings=FALSE, recursive=TRUE)


# Load library SPARQL to facilitate communication
library(SPARQL)

# The different types of queries written for R and the supporting functions
# need to be loaded:
source("rQueryFunctions.R")
source("supportingFunctions.R")


#################################################################################
# Main body starts here. This for-loop goes through each phenotype one at a time.
for(q in IterStart:IterEnd){
  
# Get starting time to keep track of each phenotype's running time
start <- Sys.time()

# Define a phenotype to test here
# phenotype <- "CD4+"
# phenotype <- "CD8+"
# phenotype <- "CD4+CD8+"
# phenotype <- "CD3 epsilon-"
 phenotype <- listPhenotypes[[q]]
# phenotype <- "CD3-CD19+CD20+IgD+CD27-"
# phenotype <- "CD45RA+CCR7+CD8+CD3+"
# phenotype <- "CD25+CCR7-"
    
cat("\n")
cat("The phenotype of interest is", phenotype, "\n")

# The following breaks up the phenotype into single markers and sorts them
# by their expression (only positive or negative expression implemented for now)
# see supportingFunctions.R for details on how 'phenoParse' works.
marker.list <- phenoParse(phenotype)

#NumberOfMarkers[[q]]      <- as.numeric(length(unlist(marker.list)))
    
# Creates another copy which will have the + and - signs. Used by flowChart and when searching files.
marker.list2 <- marker.list
# Change the . in HLA.DR to a - since the + and - signs are reserved for spilting the string
marker.list <- changeHLADR(marker.list)

# listPhenotypeOriginal[[q]]  <- phenoUnparse(marker.list)

# Put the + and- signs back on marker.list2. Used by flowChart and when searching files.
if (length(marker.list2[["Positive"]])!=0){
  for (q19 in 1:(length(marker.list2[["Positive"]]))){
    marker.list2[["Positive"]][q19] <- paste(marker.list2[["Positive"]][q19], "+", sep="")
  }}
if (length(marker.list2[["Negative"]]!=0)){
  for (q19 in 1:(length(marker.list2[["Negative"]]))){
    marker.list2[["Negative"]][q19] <- paste(marker.list2[["Negative"]][q19], "-", sep="")
  }}

# Update the marker list with the full label names from the Ontology
marker.list <- OntologyLabel(marker.list,marker.list2,markers_NickName_OntologyName)
# Make a list of the ontology names for each phenotype searched for which will be exported to a table in .csv form later
listPhenotypeUpdate[[q]]   <- phenoUnparse (marker.list)
    

skipQuery <- TRUE
res <- NULL
# Cycle through positive and negative markers.
# For positive, use a SPARQL query using the property 'has plasma membrane part'
# of each marker. Then for the negative markers, query 'lacks plasma membrane part'.
for(i in unlist(marker.list2)){
    fname <- paste(save.dirResults, 'results_', i, '.csv', sep = "")
    if(file.exists(fname)){
        tempCSV <- read.csv(fname, as.is=TRUE)
        tempCSV <- tempCSV[-c(ncol(tempCSV),ncol(tempCSV)-1)]
        res <- rbind(res, tempCSV)
    }
    else{
        cat("At least one marker was not previously queried. Querying all.\n")
        skipQuery <- FALSE
        break;
    }               
}

# If all results files exist then no querying needs to be done
if(skipQuery == TRUE){

clean.res <- tabulateResults(res)

}
else{
# Initialize result collector:
res <- NULL
for (m in marker.list[["Positive"]]){ # For each positive marker:
  cat("Locating marker", m, "\n")
  # Get relevant information about the marker -- the population names which
  # have "plasma membrane part" of this marker:
  cur.res <- queryMarker(marker = m, query.file = 'Query/hasPlasmaMembranePart.txt')
  if (nrow(cur.res) == 0){
    cat("No cell populations found which have plasma membrane part", m, "!\n")
  }

  # VERY INEFFICIENTLY, cycle through the other markers and assign a penalty 
  # according to contradiction count:
  penalties <- rep(0, nrow(cur.res))
  # Skips the following code if TRUE
  if (penalty.skip==FALSE){ # skips the penalty calculator if TRUE (very slow part of the code)
  for (other.marker in marker.list[["Negative"]]){
    for (i in 1:nrow(cur.res)){
      temp <- queryMarker(marker = other.marker, query.file = 'Query/hasPMPsingle.txt',
                          celllabel = cur.res[i, 'celllabel'])
      if (nrow(temp) > 0){
        # If one of the negatively expressed markers was found for the current
        # cell type as "has plasma membrane part", that is a contradiction which
        # should be penalized, proportionatly to the number of markers in total:
        penalties[i] <- penalties[i] - 1/length(unlist(marker.list))
      }
    }
  }
  # Similarly for the positive markers showing as negative:
  for (other.marker in setdiff(marker.list[["Positive"]], m)){
    for (i in 1:nrow(cur.res)){
      temp <- queryMarker(marker = other.marker, query.file = 'Query/lacksPMPsingle.txt',
                          celllabel = cur.res[i, 'celllabel'])
      if (nrow(temp) > 0){
        # If one of the other positively expressed markers was found for the current
        # cell type as "lacks plasma membrane part", that is a contradiction which
        # should be penalized, proportionatly to the number of markers in total:
        penalties[i] <- penalties[i] - 1/length(unlist(marker.list))
      }
    }
  }
  } # end of penalty skip if statement
  cur.res <- cbind(cur.res, penalties)
  temp.clean.res <- tabulateResults(cur.res)
  temp.name <- marker.list2[["Positive"]][which(m==marker.list[["Positive"]])]
  fname <- paste(save.dirResults, 'results_', temp.name, '.csv', sep = "")
  write.csv(temp.clean.res, fname, row.names=FALSE)
  res <- rbind(res, cur.res)
}

# Now cycle through the negative markers and construct the information similarly:
for (m in marker.list[["Negative"]]){
  cat("Locating marker", m, "\n")
  cur.res <- queryMarker(marker = m, query.file = 'Query/lacksPlasmaMembranePart.txt')
  if (nrow(cur.res) == 0){
    cat("No cell populations found which lack plasma membrane part", m, "!\n")
  }
  # VERY INEFFICIENTLY, cycle through the other markers and assign a penalty 
  # according to contradiction count:
  penalties <- rep(0, nrow(cur.res))
  # Skips the following code if TRUE
  if (penalty.skip==FALSE){ # skips the penalty calculator if TRUE (very slow part of the code)
  for (other.marker in setdiff(marker.list[["Negative"]], m)){
    for (i in 1:nrow(cur.res)){
      temp <- queryMarker(marker = other.marker, query.file = 'Query/hasPMPsingle.txt',
                          celllabel = cur.res[i, 'celllabel'])
      if (nrow(temp) > 0){
        # If one of the negatively expressed markers was found for the current
        # cell type as "has plasma membrane part", that is a contradiction which
        # should be penalized, proportionatly to the number of markers in total:
        penalties[i] <- penalties[i] - 1/length(unlist(marker.list))
      }
    }
  }
  for (other.marker in marker.list[["Positive"]]){
    for (i in 1:nrow(cur.res)){
      temp <- queryMarker(marker = other.marker, query.file = 'Query/lacksPMPsingle.txt',
                          celllabel = cur.res[i, 'celllabel'])
      if (nrow(temp) > 0){
        # If one of the other positively expressed markers was found for the current
        # cell type as "lacks plasma membrane part", that is a contradiction which
        # should be penalized, proportionatly to the number of markers in total:
        penalties[i] <- penalties[i] - 1/length(unlist(marker.list))
      }
    }
  }
  } # end of penalty skip if statement
  cur.res <- cbind(cur.res, penalties)
  temp.clean.res <- tabulateResults(cur.res)
  temp.name <- marker.list2[["Negative"]][which(m==marker.list[["Negative"]])]
  fname <- paste(save.dirResults, 'results_', temp.name, '.csv', sep = "")
  write.csv(temp.clean.res, fname, row.names=FALSE)
  res <- rbind(res, cur.res)
}

# For each returned owl object ID, tabulate how many times the result was 
# returned. This is essentially the "hits" part of the score -- telling us
# how many of the markers we have were matched to each population. From this,
# the penalty tally will be subtracted to get a final overall score for each
# population label.
clean.res <- tabulateResults(res)
}

# # Save a parital version of the results. Used when marker is queried again.
# # Will skip querying if this file exists.
# fname <- paste(save.dirPartial, 'resultsPartial_', phenotype, '.csv', sep = "")
# write.csv(res, fname, row.names=FALSE)
# cat('Initial query results saved in', fname, "\n")
    
# Save a full version of the results. (Same as partial plus penalties/number of hits/score)
fname <- paste(save.dirResults, 'results_', phenotype, '.csv', sep = "")
write.csv(clean.res, fname, row.names=FALSE)
cat('Initial query results saved in', fname, "\n")

# Now to identify the best match for this phenotype, use the tabulated
# phenotype match results:
# 1: Identify all parents of the matches:
# (NOTE that a lot of the queries used here could be combined to improve
# computational efficiency, Radina simply did not have the time to do that...)
parent.res <- matrix(nrow = 0, ncol = 5)
colnames(parent.res) <- c("x", "celllabel", "parent", "parentlabel", "score")

# If there were no hits the code will skip to the next interation instead of 
# proceeding and getting an error
if (nrow(clean.res) == 0){
    cat("Time elapsed:", TimeOutput(start), "\n")
    cat("Iterations at", q , "out of",length(listPhenotypes), "\n")
    next
}

# Added in to make the code faster when the person compiling only wants to know if the markers 
# are in the ontology or not and does not want the tree diagrams/flow charts
if (VisualizationSkip == FALSE){  

# clean.res <- clean.res.backup
    
# 2: Refine list of parents by focusing on highest scores:
scores <- as.numeric(as.character(clean.res[, "Score"]))
       
# If at least 3 perfect scores are found, only report them:
if (length(which(scores == length(unlist(marker.list)))) >= 3){
  cutoff.score <- length(unlist(marker.list))
  cat("Perfect matches found.\n")  
} else {
# Typically there will be less than perfect scores (due to marker missing/
# marker not declared in the population/non-typical phenotype)
  
  #cutoff.score <- getScoreCutoff(scores) #quantile(scores, 0.85)
  if (clean.res[1,7]>=4){
  cutoff.score <- as.integer(clean.res[1,7])-1
  }
  else{
    cutoff.score <- as.integer(clean.res[1,7])
  }
  cat("Using score cutoff of", cutoff.score, "\n")
}
    
backup.clean.res <- clean.res
    
# Extract highly scored parents only:
clean.res <- clean.res[which(scores >= cutoff.score), ]
    
    
if (length(which(scores >= cutoff.score))==1){
clean.res <- as.matrix(clean.res)
clean.res <- t(clean.res)
}
    
parent.res <- matrix(nrow = 0, ncol = 5)
colnames(parent.res) <- c("x", "celllabel", "parent", "parentlabel", "score")
for (i in 1:nrow(clean.res)){
    
    fname <- paste(save.dirParentsQuery, clean.res[i, "celllabel"], '.csv', sep = "")
    if (file.exists(fname)){  
    res <- read.csv( fname,as.is=TRUE)
    }
    else{ 
    res <- parentQuery(child.label = clean.res[i, "celllabel"], query.file = "Query/getParentClasses.txt")
    write.csv(res, fname, row.names=FALSE)      
    }
    
#   cat(i, "in", nrow(clean.res),"\n")
  res <- cbind(res, rep(clean.res[i, "Score"], nrow(res)))
  colnames(res)[ncol(res)] <- "Score"
  parent.res <- rbind(parent.res, res)
    
    
}
fname <- paste(save.dirParents, 'parent.res', phenotype, '.csv', sep = "")
write.csv(parent.res, fname, row.names=FALSE)
cat('Parent information saved in', fname, "\n")



summary <- table(parent.res[, 'parentlabel'])


# Uncomment the below part to reduce the number of nodes in the flowChart. 

# If there are too many overly specific (not a common parent to any other
# population) hits, and enough reliably general ones (parent to more than self),
# focus on the reliable ones:
# if (length(which(summary == 1)) > 10 && length(which(summary > 1)) > 3){
#   summary <- summary[-which(summary == 1)]
# }


scores2 <- sapply(names(summary), function(p) mean(as.numeric(as.character(parent.res[which(parent.res[, 'parentlabel'] == p), 'Score']))))

# The higher up in the tree -- the more times a population is called a parent to 
# others -- the more certain we are that the label applies (e.g. 'cell' is usually
# a parent to most population labels under investigation, and we are pretty sure
# whatever phenotype we are working with is at least a cell!)
# The higher the Score is, the more markers in our phenotype matched. Combining 
# these two measures gives an overall estimate of how specific and how reliable
# the hits are.
scored.summary <- summary/max(summary)*scores2/max(scores2)
s <- order(scores2, scored.summary, decreasing=TRUE)
#s <- sort(scored.summary, decreasing=TRUE, index.return=TRUE)$ix
  
# (Old version of visualization)    
# # First visualization, make sure we are getting expected (albeit not well organized) cell labels:
# colour <- 1
# 
# # The population names are printed in a size proportional to how commonly they were called a parent by another population
# # (also printed as the first value next to the cell label), and accompanied by the score (combination of how many markers got 
# # this population as a hit, minus the penalty for contradicting marker expression, as the second value):
# pdf(paste(save.dir, phenotype, ".pdf", sep = ""), width = 10, height = 5+10*(log10(length(s))))
# plot(1, 1, col="white", xlim = c(0, 2), ylim = c(0, 2), xlab="", ylab="", xaxt="none", yaxt="none", main = phenotype)
# for (i in 1:length(summary)){
#   text(x = 1, y = 1.9-2*(i-1)/length(summary), paste(names(summary)[s][i], " (", summary[s][i], ", ", round(scores[s][i], 2), ")", sep=""), cex=2+log10(summary[s][i]/max(summary)))
# }
# dev.off()
# 

  
    lengthOfs <- length(s)
# Create a list of the populations and their parents for visualization purposes  
# parent.analysis <- lapply(names(summary)[s], function(x) {
#   fname <- paste(save.dirParentsSingle, names(summary)[s], '.csv', sep = "")
#   res <- parentQuery(child.label = x, query.file = "Query/getParentClasses.txt")
#   cat("Query of",lengthOfs,"\n")    
#   return (res[which(is.element(res[, 'parentlabel'], names(summary))), 'parentlabel'])
# })
    parent.analysis <- NULL
    for (q1 in 1:length(s)){
         fname <- paste(save.dirParentsQuery, names(summary)[s[q1]], '.csv', sep = "")
        if (file.exists(fname)){  
        res <- read.csv( fname,as.is=TRUE)
        }
        else{ 
        res <- parentQuery(child.label = names(summary)[s[q1]], query.file = "Query/getParentClasses.txt")
        write.csv(res, fname, row.names=FALSE)      
        }

#         cat("Query",q1,"of",lengthOfs,"\n") 
        parent.analysis[length(parent.analysis)+1] <- list(res[which(is.element(res[, 'parentlabel'], names(summary))), 'parentlabel'])
    }
        
names(parent.analysis) <- names(summary)[s]
 
    
# Also create a list of the populations and their children to aid in the tree
# structure rendering
child.analysis <- lapply(names(summary)[s], function(x) { 
  children <- c()
  for (y in setdiff(names(summary), x)){
    if (is.element(x, parent.analysis[[y]])){
      children <- c(children, y)
    }
  }
  return (children)
})
names(child.analysis) <- names(summary)[s]
           

# Create and save a pdf file of a flow diagram
flowChart(child.analysis,clean.res, phenotype)

# (Old version of visualization)
# # The following can possibly be replaced with a better more accurate and readable
# # graphviz display of the tree structure:
# pdf(paste(save.dir, "tree_", phenotype, ".pdf", sep = ""), width = 10, height = max(10, length(parent.analysis)/2))
# # See supportingFunctions.R for treePlot(.) details
# coords <- treePlot(parent.analysis, child.analysis, scores)
# dev.off()

}# end of visualization if statement

# Check to see if all markers were hits.
  if (max(clean.res[,7])==(length(unlist(marker.list)))){
    listPhenotypeSuccess[[q]] <- ("Yes")
  }

# Compile the new row of the .csv file by updating all the list functions. 
r <- updateLists(clean.res)
listMarkerLabels[[q]] <-r[1]
listCellLabels[[q]]   <-r[2]
listPhenotypeID[[q]]  <-r[3]
listCellID[[q]]       <-r[4]
    
print(paste("Analysis for phenotype", phenotype, "complete."))
# Time for one marker/phenotype checked
cat("Time elapsed:", TimeOutput(start), "\n")
#TimeForMarkers[[q]]       <- as.numeric(Sys.time()) - as.numeric(start) 

# Shown so the person compiling knows how far the code has run
cat("Iterations at", q , "out of",length(listPhenotypes), "\n")

}
#################################################################################

# Creating a .csv file with the original phenotypes with a "Yes" or "No" indicating if it is in the ontology or not
# Plus a list of the markers and the cell types with their ontology IDs for the cases with the maximum number of hits

listPhenotypes <- cbind(listPhenotypes,listPhenotypeUpdate,listPhenotypeSuccess,listPhenotypeID,listMarkerLabels,listCellID,listCellLabels)
#listPhenotypes <- cbind(listPhenotypes,listPhenotypeSuccess)
fname <- paste(save.dir,'listPhenotypes.csv', sep = "")
write.csv(listPhenotypes, fname, row.names=FALSE)
cat("\n")  
cat("Total time was: ", TimeOutput(initialTime), "\n")

Rprof(NULL); summaryRprof(tmpf)

