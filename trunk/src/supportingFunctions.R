###################################################
# Using R to query SPARQL: supporting functions
# Author: Radina Droumeva (radina.droumeva@gmail.com)
# Date last modified: July 19, 2013
###################################################


###################################################
### Break up a phenotype into individual markers and 'signs'
#
# Args:
#   phenotype: a string to decompose into individual markers, e.g. "CD19+CD20-"
# Value:
#   A list of vectors: one of positively expressed markers and one of 
#   negatively expressed ones.
# Details:
#   TO DO: After discussion this should be changed to also support expression
#    levels other than 'positive' and 'negative' -- such as 'dim', and ??

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
  return (res)
}



###################################################
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
  sort.scores <- sort(as.numeric(clean.res[, "Score"]), 
                      decreasing = TRUE, index.return = TRUE)$ix
  return (clean.res[sort.scores, ])
}


###################################################
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


###################################################
# Locate peaks in a density function
# Input:
#   x should be a vector of values, typically arising from density(data)$y, or a 
# smoothed version of this
#   span should be a value between 0 and 1, corresponding to the proportion of 
# the range of values which the peak must cover in order to be considered valid.
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


####################################################
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


####################################################
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
