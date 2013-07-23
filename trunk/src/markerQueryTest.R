##############################################################
# Runner script for flowCL which tests functionality -- this would be the 
# script to convert into an R package eventually
# Author: Radina Droumeva (radina.droumeva@gmail.com)
# Date last modified: July 19, 2013
#
# The idea is to start with a phenotype, such as "CD4+" or "CD4+CD8+", and get 
# back the relevant portions of the ontology tree.
# There is a preliminary scoring system in place which should be assessed and 
# changed appropriately to aid in selecting the best cell label match.

# Define directory for storing results:
save.dir <- 'results/'
dir.create(save.dir, showWarnings=FALSE, recursive=TRUE)

# Get starting time to keep track of running time
start <- Sys.time()

# Define a phenotype to test here
phenotype <- "CD3-"
# phenotype <- "CD4+CD8+"
# phenotype <- "CD3 epsilon-"

# Load library SPARQL to facilitate communication
library(SPARQL)

# The different types of queries written for R and the supporting functions
# need to be loaded:
source("rQueryFunctions.R")
source("supportingFunctions.R")

# The following breaks up the phenotype into single markers and sorts them
# by their expression (only positive or negative expression implemented for now)
# see supportingFunctions.R for details on how 'phenoParse' works.
marker.list <- phenoParse(phenotype)
cat("The phenotype of interest is", phenotype, "\n")

# Cycle through positive and negative markers.
# For positive, use a SPARQL query using the property 'has plasma membrane part'
# of each marker. Then for the negative markers, query 'lacks plasma membrane
# part'.

# Initialize result collector:
res <- NULL
for (m in marker.list[["Positive"]]){ # For each positive marker:
  cat("\nLocating marker", m, "\n")
  # Get relevant information about the marker -- the population names which
  # have "plasma membrane part" of this marker:
  cur.res <- queryMarker(marker = m, query.file = 'hasPlasmaMembranePart.txt')
  if (nrow(cur.res) == 0){
    cat("No cell populations found which have plasma membrane part", m, "!\n")
  }
  # VERY INEFFICIENTLY, cycle through the other markers and assign a penalty 
  # according to contradiction count:
  penalties <- rep(0, nrow(cur.res))
  for (other.marker in marker.list[["Negative"]]){
    for (i in 1:nrow(cur.res)){
      temp <- queryMarker(marker = other.marker, query.file = 'hasPMPsingle.txt',
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
      temp <- queryMarker(marker = other.marker, query.file = 'lacksPMPsingle.txt',
                          celllabel = cur.res[i, 'celllabel'])
      if (nrow(temp) > 0){
        # If one of the other positively expressed markers was found for the current
        # cell type as "lacks plasma membrane part", that is a contradiction which
        # should be penalized, proportionatly to the number of markers in total:
        penalties[i] <- penalties[i] - 1/length(unlist(marker.list))
      }
    }
  }
  cur.res <- cbind(cur.res, penalties)
  res <- rbind(res, cur.res)
}

# Now cycle through the negative markers and construct the information similarly:
for (m in marker.list[["Negative"]]){
  cat("\nLocating marker", m, "\n")
  cur.res <- queryMarker(marker = m, query.file = 'lacksPlasmaMembranePart.txt')
  if (nrow(cur.res) == 0){
    cat("No cell populations found which lack plasma membrane part", m, "!\n")
  }
  # VERY INEFFICIENTLY, cycle through the other markers and assign a penalty 
  # according to contradiction count:
  penalties <- rep(0, nrow(cur.res))
  for (other.marker in setdiff(marker.list[["Negative"]], m)){
    for (i in 1:nrow(cur.res)){
      temp <- queryMarker(marker = other.marker, query.file = 'hasPMPsingle.txt',
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
      temp <- queryMarker(marker = other.marker, query.file = 'lacksPMPsingle.txt',
                          celllabel = cur.res[i, 'celllabel'])
      if (nrow(temp) > 0){
        # If one of the other positively expressed markers was found for the current
        # cell type as "lacks plasma membrane part", that is a contradiction which
        # should be penalized, proportionatly to the number of markers in total:
        penalties[i] <- penalties[i] - 1/length(unlist(marker.list))
      }
    }
  }
  cur.res <- cbind(cur.res, penalties)
  res <- rbind(res, cur.res)
}


# For each returned owl object ID, tabulate how many times the result was 
# returned. This is essentially the "hits" part of the score -- telling us
# how many of the markers we have were matched to each population. From this,
# the penalty tally will be subtracted to get a final overall score for each
# population label.
clean.res <- tabulateResults(res)
fname <- paste(save.dir, 'results_', phenotype, '.csv', sep = "")
write.csv(clean.res, fname, row.names=FALSE)
cat('Initial query results saved in', fname, "\n")

# Now to identify the best match for this phenotype, use the tabulated
# phenotype match results:
# 1: Identify all parents of the matches:
# (NOTE that a lot of the queries used here could be combined to improve
# computational efficiency, Radina simply did not have the time to do that...)
parent.res <- matrix(nrow = 0, ncol = 5)
colnames(parent.res) <- c("x", "celllabel", "parent", "parentlabel", "score")
for (i in 1:nrow(clean.res)){
  res <- parentQuery(child.label = clean.res[i, "celllabel"], 
                      query.file = "getParentClasses.txt")
  res <- cbind(res, rep(clean.res[i, "Score"], nrow(res)))
  colnames(res)[ncol(res)] <- "Score"
  parent.res <- rbind(parent.res, res)
}
fname <- paste(save.dir, 'parent.res', phenotype, '.csv', sep = "")
write.csv(parent.res, fname, row.names=FALSE)
cat('Parent information saved in', fname, "\n")

# 2: Refine list of parents by focusing on highest scores:
scores <- as.numeric(as.character(parent.res[, "Score"]))
# If at least 3 perfect scores are found, only report them:
if (length(which(scores == length(unlist(marker.list)))) >= 3){
  cutoff.score <- length(unlist(marker.list))
  cat("Pefect matches found.\n")  
} else {
# Typically there will be less than perfect scores (due to marker missing/
# marker not declared in the population/non-typical phenotype)
  
  cutoff.score <- getScoreCutoff(scores) #quantile(scores, 0.85)
  cat("Using score cutoff of", cutoff.score, "\n")
}

# Extract highly scored parents only:
refined.parents <- parent.res[which(scores >= cutoff.score), ]

summary <- table(refined.parents[, 'parentlabel'])
# If there are too many overly specific (not a common parent to any other
# population) hits, and enough reliably general ones (parent to more than self),
# focus on the reliable ones:
if (length(which(summary == 1)) > 10 && length(which(summary > 1)) > 3){
  summary <- summary[-which(summary == 1)]
}
scores <- sapply(names(summary), function(p) mean(as.numeric(as.character(refined.parents[which(refined.parents[, 'parentlabel'] == p), 'Score']))))

# The higher up in the tree -- the more times a population is called a parent to 
# others -- the more certain we are the label applies (e.g. 'cell' is usually
# a parent to most population labels under investigation, and we are pretty sure
# whatever phenotype we are working with is at least a cell!)
# The higher the Score is, the more markers in our phenotype matched. Combining 
# these two measures gives an overall estimate of how specific and how reliable
# the hits are.
scored.summary <- summary/max(summary)*scores/max(scores)
s <- order(scores, scored.summary, decreasing=TRUE)
#s <- sort(scored.summary, decreasing=TRUE, index.return=TRUE)$ix

# First visualization, make sure we are getting expected (albeit not well organized) cell labels:
colour <- 1

# The population names are printed in a size proportional to how commonly they were called a parent by another population (also printed as the first value next to the cell label), and accompanied by the score (combination of how many markers got this population as a hit, minus the penalty for contradicting marker expression, as the second value):
pdf(paste(save.dir, phenotype, ".pdf", sep = ""), width = 10, height = 5+10*(log10(length(s))))
plot(1, 1, col="white", xlim = c(0, 2), ylim = c(0, 2), xlab="", ylab="", xaxt="none", yaxt="none", main = phenotype)
for (i in 1:length(summary)){
  text(x = 1, y = 1.9-2*(i-1)/length(summary), paste(names(summary)[s][i], " (", summary[s][i], ", ", round(scores[s][i], 2), ")", sep=""), cex=2+log10(summary[s][i]/max(summary)))
}
dev.off()

# Create a list of the populations and their parents for visualization purposes  
parent.analysis <- lapply(names(summary)[s], function(x) {
  res <- parentQuery(child.label = x, query.file = "getParentClasses.txt")
  return (res[which(is.element(res[, 'parentlabel'], names(summary))), 'parentlabel'])
})
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

# The following can possibly be replaced with a better more accurate and readable
# graphviz display of the tree structure:
pdf(paste(save.dir, "tree_", phenotype, ".pdf", sep = ""), width = 10, height = max(10, length(parent.analysis)/2))
# See supportingFunctions.R for treePlot(.) details
coords <- treePlot(parent.analysis, child.analysis, scores)
dev.off()
print(paste("Analysis for phenotype", phenotype, "complete."))
cat("Time elapsed:", Sys.time() - start, "\n")
cat('\n')

