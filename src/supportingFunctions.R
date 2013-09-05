#################################################################################
# Using R to query SPARQL: supporting functions
# Author: Radina Droumeva (radina.droumeva@gmail.com)
# Date last modified: July 19, 2013
# Edited by: Justin Meskas (jmeskas@bccrc.ca)
# Date last modified: Sep 05, 2013
#################################################################################

#################################################################################
# Change the . in HLA.DR to a - since the + and - signs are reserved for spilting the string
changeHLADR <- function(marker.list){
if (length(marker.list[["Positive"]]) >= 1){
  for (q3 in 1:length(marker.list[["Positive"]])){
    if (marker.list[["Positive"]][q3] == "HLA.DR"){
      marker.list[["Positive"]][q3] <- "HLA-DR"
    }
  }
}
if (length(marker.list[["Negative"]]) >= 1){
  for (q3 in 1:length(marker.list[["Negative"]])){
    if (marker.list[["Negative"]][q3] == "HLA.DR"){
      marker.list[["Negative"]][q3] <- "HLA-DR"
    }
  }
}
return(marker.list)
}

#################################################################################
# Convert colour name to hexadecimal, and optionally add transparency via the 
# alpha parameter. For example, 'col2hex('red', '55') will generate a transaprent
# red. Useful when placing a legend over a plot, can make the background colour 
# of the legend col2hex("white", 70) for example.
col2hex <- function(colour, alpha = "FF"){
    colour <- col2rgb(colour)
    hex <- rgb(red = colour['red', ]/255, 
             green = colour['green', ]/255, 
              blue = colour['blue', ]/255)
    hex <- paste(hex, alpha, sep = "")
    return (hex)
}

#################################################################################
# Removes all the non first generation children from the child.analysis to produce
# a flow chart.

flowChart <- function(child.analysis, clean.res, phenotype){
  
  library("Rgraphviz")
  
  # Sort the child.analysis by starting with the cell population which has
  # the most children to the one with the least-- i.e. is the most likely to be the 
  # root parent, such as 'cell' or 'native cell':
  child.lengths <- unlist(lapply(child.analysis, length))
  sort.child <- child.analysis[order(child.lengths, decreasing=TRUE)]
  labels <- names(sort.child)
  arrows.colour <- NULL
  arrows.label <- NULL
  arrows.dashed.solid <- NULL
  arrows.head <- NULL
  
  if (ShowOntologyMarkerNames==TRUE){marker.list2 <- marker.list}
#   
#   arrows.dashed.solid <- c(labels,"temp", "temp")
#   
#   for (q20 in 1:length(arrows.dashed.solid)){
#     arrows.dashed.solid[q20] <- "solid"
#   }
#   
  
  # A quadruple for loop (probably not the most elegant, but it is still fast). Start with q8,
  # which selects the label with the most children first. q9 stores an element from the q8 label. q10 
  # and q11 look through all other elements in the q8 label and detects if the stored 
  # element (q9) is a parent of other elements in the q8 label. If this is true, then the element
  # q11 is removed from the q8 label.
  # In short, something that is not a first generation child will be removed from the list.
  for (q8 in 1:length(labels)){
    deletePoints <- NULL
    for (q9 in 1:length(child.analysis[labels[q8]][[1]])){
      
      if (is.null(child.analysis[labels[q8]][[1]][q9])){next}
      
      temp.label <- child.analysis[labels[q8]][[1]][q9]
      
      for (q10 in 1:length(child.analysis[temp.label][[1]])){
        for (q11 in 1:length(child.analysis[labels[q8]][[1]])){
          
          
          if (!is.null(child.analysis[temp.label][[1]][q10])){
            if((child.analysis[temp.label][[1]][q10] == child.analysis[labels[q8]][[1]][q11]) & (q9!=q11)){
              # stores points for removal
              deletePoints <- c(deletePoints,q11)
            }
          }
        }
      }
    }
          if(!is.null(deletePoints)){
            # removes all non first generation children
            child.analysis[labels[q8]][[1]] <- child.analysis[labels[q8]][[1]][-deletePoints] 
          }
  }
        
        # Sets up Rgraphviz with all the node titles in "labels"
        rEG <- new("graphNEL", nodes=c(labels, marker.list2[["Positive"]], marker.list2[["Negative"]]), edgemode="directed")
        
        for(q13 in 1:length(labels)){
          for (q14 in 1:length(child.analysis[labels[q13]][[1]])){
            if (!is.null(child.analysis[labels[q13]][[1]][q14])){
              # Adds an arrow from parent to child
              rEG <- addEdge(child.analysis[labels[q13]][[1]][q14], labels[q13], rEG, 1)
              arrows.colour <- c(arrows.colour, "black")
              arrows.label  <- c(arrows.label , paste(child.analysis[labels[q13]][[1]][q14], "~", labels[q13], sep=""))
              arrows.dashed.solid <- c(arrows.dashed.solid, "solid")    
              arrows.head <- c(arrows.head, "open")     
              
              
              
            }
          }
        }
        
  labels.colour <- labels
  
  # Load all phenotypes
  ColoursList <- read.csv("listOfColours.csv", header=F)
  ColoursList <- as.matrix(ColoursList)
  
  count.markers <- 0
  for (q15 in 1:length(labels)){
    for (q16 in 1:length(clean.res[,2])){
      if(clean.res[q16,2]==(labels[q15])){
        if (length(marker.list[["Positive"]])!=0){
        for (q17 in 1:length(marker.list[["Positive"]])){
          if(grepl(marker.list[["Positive"]][q17],clean.res[q16,5])==TRUE){
            count.markers <- count.markers + 1
            rEG <- addEdge(marker.list2[["Positive"]][q17],labels[q15], rEG, 1)
            arrows.colour <- c(arrows.colour, ColoursList[q17])
            arrows.label  <- c(arrows.label , paste(marker.list2[["Positive"]][q17],"~",labels[q15], sep=""))
            arrows.dashed.solid <- c(arrows.dashed.solid, "dashed")       
            arrows.head <- c(arrows.head, "none")     
            
        }}}
        if (length(marker.list[["Negative"]])!=0){
        for (q18 in 1:length(marker.list[["Negative"]])){
          if(grepl(marker.list[["Negative"]][q18],clean.res[q16,5])==TRUE){
            count.markers <- count.markers + 1
            rEG <- addEdge(marker.list2[["Negative"]][q18],labels[q15], rEG, 1)
            arrows.colour <- c(arrows.colour, ColoursList[length(marker.list2[["Positive"]])+q18])
            arrows.label  <- c(arrows.label , paste(marker.list2[["Negative"]][q18],"~",labels[q15], sep=""))
            arrows.dashed.solid <- c(arrows.dashed.solid, "dashed")
            arrows.head <- c(arrows.head, "none")     
            
        }}}
        # make perfect matches green  
        if(count.markers==length(unlist(marker.list))){        
        labels.colour[q15] <- "lightgreen"
        count.markers <- 0
        }
        else{ # colour the partial matches a beige colour
          labels.colour[q15] <- "bisque"
          count.markers <- 0     
        }
        break
      }
      else{ # colour all the non-important parents white
        labels.colour[q15] <- "white"
      }
      
    }
  }
  if (length(marker.list[["Positive"]])!=0){ # colour positive markers sky blue
  for (q19 in 1:(length(marker.list[["Positive"]]))){
  labels.colour[length(labels)+q19] <- "skyblue"
  }}
  if (length(marker.list[["Negative"]]!=0)){
  for (q19 in 1:(length(marker.list[["Negative"]]))){ #colour negative markers a light red colour
    labels.colour[length(labels)+length(marker.list[["Positive"]])+q19] <- "lightcoral"
  }}
  nAttrs <- list()
  eAttrs <- list()

#   nAttrs$color <- structure(c(labels.colour), .Names = c(labels))
  nAttrs$fillcolor <- structure(c(labels.colour), .Names = c(labels, marker.list2[["Positive"]], marker.list2[["Negative"]]))
  eAttrs$color <- structure(c(arrows.colour), .Names = c(arrows.label))
  
#  edgesInfo <- buildEdgeList(rEG, edgeAttrs=eAttrs, defAttrs=defAttrs$edge)  
  eAttrs$style <- structure(c(arrows.dashed.solid), .Names = c(arrows.label))
  eAttrs$arrowhead <- structure(c(arrows.head), .Names = c(arrows.label))
  
#   for (q21 in 1:length(arrows.label)){
#   edgesInfo[[arrows.label[q21]]]@attrs$arrowhead <- arrows.head[q21]
#   }
  
#   eAttrs$style <- c("CD4 molecule~naive regulatory T cell"="dotted")
#   edgeData(rEG, "CD4 molecule~naive regulatory T cell", c("style")) <- c("dashed")
#   edges[["CD4 molecule~naive regulatory T cell"]]@attrs$arrowhead <- "none"
#   eAttrs$arrowhead <- c("CD4 molecule~naive regulatory T cell"="none")
  
  attrs <- list( node=list(shape="ellipse", fontsize = 14, fixedsize=FALSE),graph=list(rankdir="BT"))
  #, edge=list(color="green", arrowhead="none")
#   attrs$edge$style <- structure(c(arrows.dashed.solid), .Names = c(arrows.label))
#   attrs$edge$arrowhead <- structure(c(arrows.head))
 
        # create a pdf flow chart
        child.file.name <- paste(save.dir, "tree_", phenotype, ".pdf", sep = "")     
        pdf(file=child.file.name, width=10.5, height=8)    
        plot(rEG, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs) 
  
#   plot(rEG, attrs=list(node=list(shape="ellipse", fontsize = 14, fixedsize=FALSE,label="foo", fillcolor="lightgreen"),
#                        edge=list(color="green"),
#                        graph=list(rankdir="LR")))     
        #plot(rEG, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs) 
        dev.off()
        
}

#################################################################################
# Locate peaks in a density function
# Input:
#   x should be a vector of values, typically arising from density(data)$y, or a 
# smoothed version of this
#   span should be a value between 0 and 1, corresponding to the proportion of 
# the range of values which the peak must cover in order to be coneural cell adhesion molecule 1nsidered valid.
getPeaks <- function(x, span=.05){
  x <- as.vector(x)
  peaks <- c()
  span.len <- max(ceiling(length(x)*span/2), 1)

  for (i in 2:(length(x) - 1)){
    if (x[i-1] < x[i] && x[i+1] <= x[i]){ # if a peak in its immediate neighbourhood
      close.peak <- which(is.element(seq(from=i-span.len, to=i+span.len, by=1), peaks))[1] #locate closest peak within span
      if (!is.na(close.peak)){ # There should never be more than 1, as we are sweeping from left to right and replacing close peaks with the maximum in the span-neighbourhood as we go
        if(x[close.peak] < x[i]) # if current peak is bigger than the pre-existing close peak, replace it
          peaks[which(peaks == close.peak)] <- i
      } else {
        if (!any(x[max(i-span.len, 1):(i-1)] > x[i]) && !any(x[(i+1):min(i+span.len, length(x))] > x[i]) && (i > span.len)){
          peaks <- c(peaks, i)
        }
      }
    }
  }
  res <- rep(FALSE, length(x))
  res[peaks] <- TRUE
  return (res)
}


#################################################################################
# Whenever there is no perfect score, instead of trying to use the quantile to set
# a threshold, use the density of scores and grab the highest peak:
getScoreCutoff <- function(scores){
  dens <- density(scores)
  peaks <- which(getPeaks(dens$y))
  peaks <- peaks[which(dens$y[peaks] > 0.05*max(dens$y))]
  if (length(peaks) <= 1)
    return (dens$x[1])
  #
  peaks <- tail(peaks, 2)
  min.point <- which.min(dens$y[peaks[1]:peaks[2]])
  cutoff <- dens$x[peaks[1]:peaks[2]][min.point]
  
  # To visualize, run the following:
  # plot(dens)
  # abline(v = cutoff)
  
  return (cutoff)
}

#################################################################################
### Break up a phenotype into individual markers and 'signs'
#
# Args:
#   phenotype: a string to decompose into individual markers, e.g. "CD19+CD20-"
# Value:
#   A list of vectors: one of positively expressed markers and one of 
#   negatively expressed ones.
# Details:
#   TO DO: After discussion this should be changed to also support expression
#    levels other than 'positive' and 'negative' -- such as 'dim', 'low', 'lo', 'high' and ??

phenoParse <- function(phenotype) {
  if (!is.character(phenotype)){
    warning("Phenotype is not a valid string!")
    try (phenotype <- as.character(phenotype))
  }
  
  # First split up the string based on + or - to get the markers
  markers <- unlist(strsplit(x = phenotype, split = "\\+|\\-"))
  # Next, split the original string based on the markers found above, leaving
  # only the signs (remove first result as there is no leading sign)
  
  signs <- unlist(strsplit(x = phenotype, 
                           split = paste(markers, sep = "", collapse="|")))[-1]
  
  # Return a list of positive and negative markers
  res <- list(`Positive` = markers[signs == "+"],
              `Negative` = markers[signs == "-"])
  # if (strpos(phenotype,'HLA-DR'){
  
  #  }
  return (res)
}



#################################################################################
# Creates a string with the updated phenotypes to have better formatting for the .csv file
phenoUnparse <- function(marker.list){
  temp.string <- NULL
  if(length(marker.list[["Positive"]])>=1){
    for(q6 in 1:length(marker.list[["Positive"]])){
      temp.string <-  paste(temp.string, marker.list[["Positive"]][q6],"\n",sep="")
    }}
  if(length(marker.list[["Negative"]])>=1){
    for(q6 in 1:length(marker.list[["Negative"]])){
        temp.string <-  paste(temp.string, marker.list[["Negative"]][q6],"\n",sep="")     
    }}
  # removes the \n from the last line for formatting reasons in the .csv file
  temp.string <- substr(temp.string,1,nchar(temp.string)-1) 
  return(temp.string)
}

#################################################################################
### Condense a results table to the unique entries only and tabulate repeated
### hits
#
# Args:
#   res: a matrix containing SPARQL query results
# Value:
#   A cleanly tabulated results matrix (table)
# Details:
#   TO DO: For now, this relies on the first column having the unique IDs of
#   the owl objects returned.
tabulateResults <- function(res){
  number.of.hits <- table(res[, 1])
  res <- cbind(res, number.of.hits[res[, 1]])
  colnames(res)[ncol(res)] <- "Number Of Hits"
  unique.ids <- unique(res[, 1])
  clean.res <- matrix("", nrow = length(unique.ids), ncol = ncol(res)+1)
  colnames(clean.res) <- c(colnames(res), "Score")
  rownames(clean.res) <- unique.ids
  for (id in unique.ids){
    locate.entries <- which(res[, 1] == id)
    clean.res[id, ] <- c(id, res[locate.entries[1], 2], 
                         apply(res[locate.entries, 3:(ncol(res) - 2)], 2, paste, 
                               collapse = "\n"), sum(res[locate.entries, 'penalties']),
                         res[locate.entries[1], ncol(res)], 0)
    clean.res[id, "Score"] <- as.numeric(clean.res[id, "Number Of Hits"]) + 
      as.numeric(clean.res[id, "penalties"])
  }
  # Required because of a bug cause from only having 1 row (f[1,] notation was a problem)
  if (nrow(clean.res)>=2){
    sort.scores <- sort(as.numeric(clean.res[, "Score"]), 
                        decreasing = TRUE, index.return = TRUE)$ix
    return (clean.res[sort.scores, ])
  }
  else{
    return (clean.res)
  }
}

#################################################################################
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
TimeOutput(Sys.Date())


#################################################################################
# Plot results in a tree structure: THIS SHOULD PROBABLY BE REPLACED WITH GRAPHVIZ OR ANYTHING BETTER!!!!
treePlot <- function(parent.analysis, child.analysis, scores, xlim=NULL, ylim=NULL){
  parent.analysis <- lapply(parent.analysis, unique)
  child.analysis <- lapply(child.analysis, unique)
  
  # Start creating the tree structure at the top, by starting with the cell 
  # population which has the most children -- i.e. is the most likely to be the 
  # root parent, such as 'cell' or 'native cell':
  child.lengths <- unlist(lapply(child.analysis, length))
  parent.lengths <- unlist(lapply(parent.analysis, length))
  sort.child <- child.analysis[order(child.lengths, decreasing=TRUE)]
  labels <- names(sort.child)

  # Number of clusters indicates number of independent beginnings (roots) in the tree
  num.clusters <- max(length(which(parent.lengths == min(parent.lengths))), 1)
  cat("Number of population clusters:", num.clusters, "\n")
  # This is a rough measure of the depth of the tree, for plotting purposes:
  m <- max(child.lengths)
   # Initialize some parameters to help arrange the text in the plot well:
  colour <-1 # reset to first line to be black
  coords <- matrix(0, ncol=2, nrow=length(labels)) #keep track of where each population is plotted
  rownames(coords) <- labels
 
  cluster <- 1 # start with first cluster of relations (first root)
  signs <- rep(1, length(labels)) # alternate direction in which branching occurs
  signs[seq(2, length(signs), by=2)] <- -1
  names(signs) <- labels
  levels <- rep(0, length(labels)) # keep track of how many branches there are at the current level
  names(levels) <- labels

  # For connecting lines that need to be plotted for populations not yet plotted (i.e. cross-cluster lines)
  plot.second.round <- list()
  if (is.null(xlim)){
    xlim = c(-1/2*length(child.lengths), 2*length(child.lengths))
  }
  if (is.null(ylim)){
    ylim = c(-m/2, m)
  }
  plot(1,1, col="white", xlim = xlim, ylim = ylim)
  correction.y <- 1.1*strheight("CD") # typical text height
  x.threshold <- 1.2*strwidth("native cell") # typical text length
  cat("Correction of y: ")
  print(correction.y)
  i <- strheight("A") # default horizontal spacing between tree branches
  # Initial placement of first root is in the middle (or 1/3rd of the way if more roots) at the top
  x <- length(child.analysis)/(1+num.clusters)
  y <- m
  
  # Continue plotting until we have plotted all labels. As we plot them, we will remove them from the list
  # of still-to-be-plotted.
  while (!is.null(labels) && length(labels)!=0){
    
    # Root parent to plot:
    parent <- labels[1]
    
    # If it is in the list of hits matching our phenotype, and not just a parent of a hit, colour it red!
    if (is.element(parent, clean.res[, 'celllabel'])){
          col <- "red"
      } else {
        col <- "black"
    }
    
    # Plot the parent, record its coordinates, and remove from the list of cell populations still to plot (labels):
    text(x, y, paste(parent, " (", scores[parent], ")", sep = ""), col=col)
    coords[parent, ] <- c(x, y)
    labels <- setdiff(labels, parent)

    # Now we will plot each of this root parent's children. To maximize our chance of doing this 'right',
    # we start with the child which has the lowest number of parents -- i.e. being closer to 'root'.  
    order.children <- names(sort(parent.lengths[child.analysis[[parent]]], decreasing=FALSE))
    order.children <- intersect(labels, order.children)
    while (!is.null(order.children) && length(order.children)!=0){
      child <- order.children[1]
      # Locate all non-self parents for the current child:
      this.child.parents <- setdiff(parent.analysis[[child]], child)
      
      # To find the most direct parent, choose the one which itself has the most parents (could be more than one):
      max.parent.count <- max(parent.lengths[this.child.parents])
      direct.parents <- this.child.parents[which(parent.lengths[this.child.parents] == max.parent.count)]

      # Draw a line connecting this child's direct parents if possible, otherwise set aside to draw later:
      child.printed <- FALSE
      for (d in direct.parents){
        if (all(coords[d, ] == c(0, 0))){ # If this direct parent has not yet been plotted:
          plot.second.round[[length(plot.second.round)+1]] <- c(child, d)
        } else {
          sign <- signs[d]
          signs[d] <- -1*signs[d]
          levels[d] <- levels[d] + 1
          if (!child.printed){
            x.cur <- coords[d, 1] + levels[d]*10*i*sign
            y.cur <- coords[d, 2] - 2*correction.y
            
            condition.to.change <- TRUE
            # Check if there is something else printed with coordinates too close to current,
            # and if so, try printing lower. Continue checking until you find space that is not in the way.
            while (any(condition.to.change)){
              condition.to.change <- sapply((1:nrow(coords))[-which(rownames(coords) == child)], function(j) {
                      part1 <- (abs(coords[j, 2] - y.cur) < correction.y)
                      part2 <- (abs(coords[j, 1] - x.cur) < (1/2*strwidth(rownames(coords)[j]) + 1/2*strwidth(child) + 2*strwidth("(A)")))
                      result <- (part1 && part2)
                      return (result)
                    })
              if (any(condition.to.change)){
                y.cur <- y.cur - correction.y
              }
            }
            child.printed <- TRUE
          }
          # Double check again that the final chosen printing spot is fine.
          condition.to.change <- TRUE
          while (any(condition.to.change)){
            condition.to.change <- sapply((1:nrow(coords))[-which(rownames(coords) == child)], function(j) {
                    part1 <- (abs(coords[j, 2] - y.cur) < correction.y)
                    part2 <- (abs(coords[j, 1] - x.cur) < (1/2*strwidth(rownames(coords)[j]) + 1/2*strwidth(child) + 2*strwidth("(A)")))
                    result <- (part1 && part2)
                    return (result)
                  })
            if (any(condition.to.change)){
              y.cur <- y.cur - correction.y
            }
          }
          lines(c(coords[d, 1], x.cur), c(coords[d, 2], y.cur), col = col2hex(colour, "55"), lwd=3)
          points(x.cur, y.cur, pch=19, col = col2hex(colour, "33"), cex = 2)
          colour <- colour+1   
        }
      }
      condition.to.change <- TRUE
      while (any(condition.to.change)){
        condition.to.change <- sapply((1:nrow(coords))[-which(rownames(coords) == child)], function(j) {
                part1 <- (abs(coords[j, 2] - y.cur) < correction.y)
                part2 <- (abs(coords[j, 1] - x.cur) < (1/2*strwidth(rownames(coords)[j]) + 1/2*strwidth(child) + 2*strwidth("(A)")))
                result <- (part1 && part2)
                return (result)
              })
        if (any(condition.to.change)){
          y.cur <- y.cur - correction.y
        }
      }
      # Print the child's name, finally:
      x <- x.cur
      y <- y.cur
      if (is.element(child, clean.res[, 'celllabel'])){
        col <- "red"
      } else {
        col <- "black"
      }
      text(x, y, paste(child, " (", scores[child], ")", sep = ""), col=col)
      coords[child, ] <- c(x, y)
      labels <- setdiff(labels, child)
      order.children <- intersect(order.children, labels)
    }
    # This loop ends when the first root population has all its children printed.
    cluster <- cluster + 1
    x <- length(child.analysis)*cluster/(1+num.clusters)
    y <- m
  }
  
  print(plot.second.round)
  if (length(plot.second.round) > 0){
    for (i in 1:length(plot.second.round)){
      colour <- colour+1
      pts1 <- coords[plot.second.round[[i]][1], ]
      pts2 <- coords[plot.second.round[[i]][2], ]
      lines(c(pts1[1], pts2[1]), c(pts1[2], pts2[2]), col = col2hex(colour, "55"), lwd=3)
      points(pts2, pch = 19, col=col2hex(colour, "33"), cex=2)
    }
  }
  return (coords)
} 

#################################################################################
# Finds and stores information for display in a .csv file
updateLists <- function(clean.res){
  # Creates the lists MarkerLabels, CellLabels, PhenotypeID and CellID
  # Column 7 is the Number of Hits
  # Column 5 is the type of markers that were hits
  # Column 2 is the type of cell with the corresponding markers
  # Column 4 is the marker ID
  # Column 1 is the cell ID
  # The code could be rewritten to look for column names instead of relying on thecolumn numbers
  BreakTrue <- FALSE
  if (max(clean.res[,7])>=1){
    for (q2 in 1:min(length(clean.res[,7]),5)){
      if (q2 == 1 & (clean.res[q2,7]) == max(clean.res[,7])){
        listMarkerLabels.temp <- paste(q2,") ",(clean.res[q2,5]),"\n", sep="") 
        listCellLabels.temp   <- paste(q2,") ",(clean.res[q2,2]),"\n", sep="")
        
        # Look in the Label of the marker ID for "PR" then pulls out the PR and the numbers that follow 
        temp.index <- gregexpr("PR", (clean.res[q2,4]))
        temp.string <- ""
        if ((temp.index[[1]][1]!=-1)){
          for (q7 in 1:length(temp.index[[1]])){
            temp.string <- paste(temp.string, substr(clean.res[q2,4],temp.index[[1]][q7],temp.index[[1]][q7]+11),"\n", sep="")
            if(q7 == length(temp.index[[1]])){temp.string <- substr(temp.string,1,nchar(temp.string)-1)} # Removes the last \n 
          }}
        # Look in the Label of the marker ID for "GO" then pulls out the GO and the numbers that follow
        temp.index <- gregexpr("GO", (clean.res[q2,4]))
        if ((temp.index[[1]][1]!=-1)){
          if (temp.string!=""){
            temp.string <- paste(temp.string,"\n",sep="")
          }
          for (q7 in 1:length(temp.index[[1]])){
            temp.string <- paste(temp.string, substr(clean.res[q2,4],temp.index[[1]][q7],temp.index[[1]][q7]+9),"\n", sep="")
            if(q7 == length(temp.index[[1]])){temp.string <- substr(temp.string,1,nchar(temp.string)-1)} # Removes the last \n      
          }}
        # Stores this string into the lists
        listPhenotypeID.temp  <- paste(q2,") ",temp.string,"\n", sep="")
        
        # Look in the Label of the Cell ID for "CL" then pulls out the CL and the numbers that follow
        temp.index <- gregexpr("CL", (clean.res[q2,1]))
        temp.string <- ""
        if ((temp.index[[1]][1]!=-1)){
          for (q7 in 1:length(temp.index[[1]])){
            temp.string <- paste(temp.string, substr(clean.res[q2,1],temp.index[[1]][q7],temp.index[[1]][q7]+9),"\n", sep="")
            if(q7 == length(temp.index[[1]])){temp.string <- substr(temp.string,1,nchar(temp.string)-1)} # Removes the last \n 
          }}
        # Stores this string into the lists
        listCellID.temp  <- paste(q2,") ",temp.string,"\n", sep="")
      }
      if (q2 > 1 & (clean.res[q2,7]) == max(clean.res[,7])){
        listMarkerLabels.temp <- paste(listMarkerLabels.temp , q2, ") ",(clean.res[q2,5]),"\n", sep="")
        listCellLabels.temp   <- paste(listCellLabels.temp   , q2, ") ",(clean.res[q2,2]),"\n", sep="")
        
        # Look in the Label of the marker ID for "PR" then pulls out the PR and the numbers that follow 
        temp.index <- gregexpr("PR", (clean.res[q2,4]))
        temp.string <- ""
        if ((temp.index[[1]][1]!=-1)){
          for (q7 in 1:length(temp.index[[1]])){
            temp.string <- paste(temp.string, substr(clean.res[q2,4],temp.index[[1]][q7],temp.index[[1]][q7]+11),"\n", sep="")
            if(q7 == length(temp.index[[1]])){temp.string <- substr(temp.string,1,nchar(temp.string)-1)} # Removes the last \n 
          }}
        # Look in the Label of the marker ID for "GO" then pulls out the GO and the numbers that follow
        temp.index <- gregexpr("GO", (clean.res[q2,4]))
        if ((temp.index[[1]][1]!=-1)){
          if (temp.string!=""){
            temp.string <- paste(temp.string,"\n",sep="")
          }
          for (q7 in 1:length(temp.index[[1]])){
            temp.string <- paste(temp.string, substr(clean.res[q2,4],temp.index[[1]][q7],temp.index[[1]][q7]+9),"\n", sep="")
            if(q7 == length(temp.index[[1]])){temp.string <- substr(temp.string,1,nchar(temp.string)-1)} # Removes the last \n        
          }}
        # Stores this string into the lists
        listPhenotypeID.temp  <- paste(listPhenotypeID.temp  , q2, ") ",temp.string,"\n", sep="")
        
        # Look in the Label of the Cell ID for "CL" then pulls out the CL and the numbers that follow
        temp.index <- gregexpr("CL", (clean.res[q2,1]))
        temp.string <- ""
        if ((temp.index[[1]][1]!=-1)){
          for (q7 in 1:length(temp.index[[1]])){
            temp.string <- paste(temp.string, substr(clean.res[q2,1],temp.index[[1]][q7],temp.index[[1]][q7]+9),"\n", sep="")
            if(q7 == length(temp.index[[1]])){temp.string <- substr(temp.string,1,nchar(temp.string)-1)} # Removes the last \n 
          }}
        # Stores this string into the lists
        listCellID.temp  <- paste(listCellID.temp  , q2, ") ",temp.string,"\n", sep="")
      }
      if ((clean.res[q2,7]) != max(clean.res[,7])){
        BreakTrue <- TRUE
        break
      }  
    }# end of for loop
  }# end of if statement
  
  # If there is more than 5 elements to store, then the rest is cut off and a "+ more" is displayed
  if (length(clean.res[,7])>5 & BreakTrue == FALSE){
    listMarkerLabels.temp <- paste(listMarkerLabels.temp, "+ more")
    listCellLabels.temp   <- paste(listCellLabels.temp,   "+ more")
    listPhenotypeID.temp  <- paste(listPhenotypeID.temp,  "+ more")
    listCellID.temp       <- paste(listCellID.temp,       "+ more")
  }
  # Removes the last \n from the string for better formatting int the .csv file
  else{
    listMarkerLabels.temp <- substr(listMarkerLabels.temp,1,nchar(listMarkerLabels.temp)-1)
    listCellLabels.temp   <- substr(listCellLabels.temp,  1,nchar(listCellLabels.temp)-1)
    listPhenotypeID.temp  <- substr(listPhenotypeID.temp, 1,nchar(listPhenotypeID.temp)-1)
    listCellID.temp       <- substr(listCellID.temp,      1,nchar(listCellID.temp)-1)
  }
  
  return(c(listMarkerLabels.temp, listCellLabels.temp, listPhenotypeID.temp, listCellID.temp))
}


