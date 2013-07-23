###################################################
# Using R to query SPARQL: query functions
# Author: Radina Droumeva (radina.droumeva@gmail.com)
# Last date modified: July 19/2013
###################################################


###################################################
### Get ID and label for single marker search
#
# Args:
#   marker: a string to search for, such as "CD19"
#   query.file: a .txt file name (full path if not in current workspace
#   directory) containing the SPARQL query.
#   celllabel: a string matching the label of a cell type for which we are 
#   running the query. This is applicable for certain queries only!
#   Specifically, for hasPMPsingle.txt and lacksPMPsingle.txt, which
#   check if population 'celllabel' has or lacks plasma membrane part marker 'm'
# Value:
#   res: an Nx2 matrix containing all unique pairs of (ID, label) such that
#   the passed 'marker' is matched by regular expression to the selected labels
#   Note that the IDs need not be unique, but the pairs of (ID, label) are, or
#   an Nx5 matrix for unique ID, label, parent label, marker, marker label, if 
#   'celllabel' is not null.
# Details:
#   At this stage the key regular expression is of the form:
#    "CD19$|CD19[^0-9a-zA-Z][^ alpha][^ beta]".
#   This selects entries where the match is either of the exact marker passed 
#   followed by the end of line, or followed by anything other than a number
#   or letter or space and 'alpha', or space and 'beta'. This ensures scenarios 
#   such as the following do not occur:
#   Case 1: Looking for CD19 returns CD19, CD193, etc.
#   Case 2: Looking for Ly6 returns Ly6g, etc.
#   Case 3: Looking for CD8 returns CD8 alpha chain

queryMarker <- function(marker = NULL, query.file = "getMatchingSynonyms.txt",
                        celllabel = NULL){

  # TO DO: improve input check
  # For now, a NULL, NA or length 0 phenotype will have an empty results
  #  matrix returned.
  if (!is.character(marker) || length(marker) == 0){
    warning("No marker found.")
    res <- matrix(ncol = 2, nrow = 0)
    colnames(res) <- c("ID", "Synonym Match")
    return (res)
  }

  # Define the cell.ctde.net SPARQL endpoint
  endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"

  # The options for 'read.csv' ensure that lines remain strings and not factors,
  # the lines are parsed based on the new line character (not a comma),
  # any quotes found in the parsing of the text are left intact,
  # and none of the lines being read should be treated as a header.
  q <- read.csv(query.file, as.is = TRUE, sep = "\n", quote = "", 
                header = FALSE)[, 1]
  prefix.info <- read.csv('prefixInfo.txt', as.is = TRUE, sep = "\n",
                          quote = "", header = FALSE)[, 1]
  # Concatenate the query preceeded by all prefix information as a single string
  # to be passed to SPARQL
  query <- paste(c(prefix.info, q), collapse="\n")

  # Prepare marker for query by ensuring the marker is either followed by the  
  # end of the line or it has a symbol other than a letter or number after it, 
  # as that may indicate a different marker
  marker <- paste(marker, "$|", marker, "[^0-9a-zA-Z-+][^ alpha][^ beta]", collapse="", sep="")
  query <- gsub("\\$marker", marker, query)
  if (!is.null(celllabel)){
    celllabel <- paste("^", celllabel, "$", sep = "", collapse = "")
    query <- gsub("\\$celllabel", celllabel, query)
  }

  # Execute query
  res <- SPARQL(url=endpoint, query)$results
  return (res)
}

# Similar to above, but for finding the parents for a specific cell population
# by cell label.
# NOTE: I tried matching by ID and for some reason I could do it on the HE group
# but not here for some reason, that's why I'm matching by label...
parentQuery <- function(child.label = "common myeloid progenitor", 
                         query.file = "getParentClasses.txt"){

  # Define the cell.ctde.net SPARQL endpoint
  endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"

  # The options for 'read.csv' ensure that lines remain strings and not factors,
  # the lines are parsed based on the new line character (not a comma),
  # any quotes found in the parsing of the text are left intact,
  # and none of the lines being read should be treated as a header.
  q <- read.csv(query.file, as.is = TRUE, sep = "\n", quote = "", 
                header = FALSE)[, 1]
  prefix.info <- read.csv('prefixInfo.txt', as.is = TRUE, sep = "\n",
                          quote = "", header = FALSE)[, 1]
  # Concatenate the query preceeded by all prefix information as a single string
  # to be passed to SPARQL
  query <- paste(c(prefix.info, q), collapse="\n")
  # Add "^" to beginning and "$" to end to find an exact match for the label
  child.label <- paste("^", child.label, "$", sep = "", collapse = "")
  query <- gsub("\\$label", child.label, query)
  # Execute query
  res <- SPARQL(url=endpoint, query)$results 
  return (res)
}



