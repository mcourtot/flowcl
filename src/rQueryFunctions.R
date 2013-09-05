#################################################################################
# Using R to query SPARQL: query functions
# Author: Radina Droumeva (radina.droumeva@gmail.com)
# Last date modified: July 19, 2013
# Edited by: Justin Meskas (jmeskas@bccrc.ca)
# Date last modified: Sept 05, 2013
#################################################################################


#################################################################################
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
#    "CD19$|CD19[^0-9a-zA-Z][^ alpha][^ beta][^ epsilon]".
#   This selects entries where the match is either of the exact marker passed 
#   followed by the end of line, or followed by anything other than a number
#   or letter or space and 'alpha', or space and 'beta'. This ensures scenarios 
#   such as the following do not occur:
#   Case 1: Looking for CD19 returns CD19, CD193, etc.
#   Case 2: Looking for Ly6 returns Ly6g, etc.
#   Case 3: Looking for CD8 returns CD8 alpha chain

# The above mentioned "regular expression" was incorrect. The [^ alpha][^ beta][^ epsilon] command excluded the words 
# that had one of the characters a, l, p, h and a in the first letter after the space, the character b,e,t and a in the second
# letter after the space and the character e, p, s, i, l, o and n in the third letter. The updated command [^ abe]
# will exclude the words starting with an a, b or e. This will exclude alpha, beta and epsilon, however, it will exclude any 
# other words starting with a, b or e. (A potential bug waiting to happen) Research into proper use of Regular Expression" yielded
# a \b command (i.e. \bbeta\b). Unfortunately this did not work.

#################################################################################
# function used by OntologyLabel
applyOntologyExceptionsNegative1 <- function(marker.list,q1){ 
    
  updatedException <- FALSE
  if(marker.list[q1]=="CD8")   {updatedException <- TRUE}
  if(marker.list[q1]=="IgD")   {updatedException <- TRUE}
  if(marker.list[q1]=="CD3")   {updatedException <- TRUE}
  if(marker.list[q1]=="HLA-DR"){updatedException <- TRUE}  
  return(updatedException)
}
#################################################################################
# function used by OntologyLabel
applyOntologyExceptionsPositive1 <- function(marker.list,q1){ 
  
  updatedException <- FALSE
  if(marker.list[q1]=="CD8")   {updatedException <- TRUE}
  if(marker.list[q1]=="IgD")   {updatedException <- TRUE}
  if(marker.list[q1]=="CD3")   {updatedException <- TRUE}
  if(marker.list[q1]=="HLA-DR"){updatedException <- TRUE}  
  return(updatedException)
}
#################################################################################
# function used by OntologyLabel
applyOntologyExceptionsNegative2 <- function(marker.list,q1){ 
  
  if(marker.list[q1]=="CD8"){
    cat("Marker", marker.list[q1], "has been FORCED to updated to T cell receptor co-receptor CD8\n")
    marker.list[q1] <-"T cell receptor co-receptor CD8"
  }
  if(marker.list[q1]=="IgD"){
    cat("Marker", marker.list[q1], "has been FORCED to updated to IgD immunoglobulin complex\n")
    marker.list[q1] <-"IgD immunoglobulin complex"
  }
  if(marker.list[q1]=="CD3"){
    cat("Marker", marker.list[q1], "has been FORCED to updated to CD3 epsilon\n")
    marker.list[q1] <-"CD3 epsilon"
  }
  if(marker.list[q1]=="HLA-DR"){
    cat("Marker", marker.list[q1], "has been FORCED to updated to MHC class II protein complex\n")
    marker.list[q1] <-"MHC class II protein complex"
  }
  return(marker.list)
}
#################################################################################
# function used by OntologyLabel
applyOntologyExceptionsPositive2 <- function(marker.list,q1){ 
  
  if(marker.list[q1]=="CD8"){
    cat("Marker", marker.list[q1], "has been FORCED to updated to T cell receptor co-receptor CD8\n")
    marker.list[q1] <-"T cell receptor co-receptor CD8"
  }
  if(marker.list[q1]=="IgD"){
    cat("Marker", marker.list[q1], "has been FORCED to updated to IgD immunoglobulin complex\n")
    marker.list[q1] <-"IgD immunoglobulin complex"
  }
  if(marker.list[q1]=="CD3"){
    cat("Marker", marker.list[q1], "has been FORCED to updated to alpha-beta T cell receptor complex\n")
    marker.list[q1] <-"alpha-beta T cell receptor complex"
  }
  if(marker.list[q1]=="HLA-DR"){
    cat("Marker", marker.list[q1], "has been FORCED to updated to MHC class II protein complex\n")
    marker.list[q1] <-"MHC class II protein complex"
  }
  return(marker.list)
}
#################################################################################
# Searches the Ontology to find the correct label for each marker
OntologyLabel <- function(marker.list, marker.list2, markers_NickName_OntologyName){
  
  if(length(marker.list[["Positive"]])!=0){            
  for(q1 in 1:length(marker.list[["Positive"]])){
    
    # skips query if there is a file named "markers_NickName_OntologyName" with all the marker ontology labels
    fname <- "CL_Markers/markers_NickName_OntologyName.csv"
    if(file.exists(fname)){
    markers_NickName_OntologyName <- read.csv(fname , header=F)    
    markers_NickName_OntologyName <- as.matrix(markers_NickName_OntologyName)

    temp.location <- which(markers_NickName_OntologyName[,1]==marker.list2[["Positive"]][q1])
    if(length(temp.location)==1){
        cat("Marker", marker.list[["Positive"]][q1], "is called", markers_NickName_OntologyName[temp.location,2], "\n")
        marker.list[["Positive"]][q1] <- markers_NickName_OntologyName[temp.location,2]
        next
    }
    }
      
    # Change the short name markers in marker.list to the marker labels in the Ontology
    # This is only needed for the ones that the code has trouble finding
    # Hopefully with an updated Ontology these next three line can be deleted
    updatedException <- applyOntologyExceptionsPositive1(marker.list[["Positive"]], q1)
    marker.list[["Positive"]] <- applyOntologyExceptionsPositive2(marker.list[["Positive"]], q1)
    if ( updatedException == TRUE){next}
    
    temp.marker <- marker.list[["Positive"]][q1]
    # Define the cell.ctde.net SPARQL endpoint
    endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"
    
    # The options for 'read.csv' ensure that lines remain strings and not factors,
    # the lines are parsed based on the new line character (not a comma),
    # any quotes found in the parsing of the text are left intact,
    # and none of the lines being read should be treated as a header.
    que <- read.csv("Query/hasProperLabel.txt", as.is = TRUE, sep = "\n", quote = "", 
                  header = FALSE)[, 1]
    prefix.info <- read.csv('Query/prefixInfo.txt', as.is = TRUE, sep = "\n",
                            quote = "", header = FALSE)[, 1]
    # Concatenate the query preceeded by all prefix information as a single string
    # to be passed to SPARQL
    query <- paste(c(prefix.info, que), collapse="\n")
    
    # Prepare marker for query by ensuring the marker is either followed by the  
    # end of the line or it has a symbol other than a letter or number after it, 
    # as that may indicate a different marker
    temp.marker <- paste(temp.marker, "$|", temp.marker, "[^0-9a-zA-Z-+][^ abes]", collapse="", sep="")
    query <- gsub("\\$marker", temp.marker, query)
    # Execute query  
    res <- SPARQL(url=endpoint, query)$results
    if (nrow(res)==1){
      temp.marker <- res[1,2]          
      cat("Marker", marker.list[["Positive"]][q1], "has been updated to", temp.marker, "\n")
        
      # Update the markers_NickName_OntologyName.csv file
      fname <- "CL_Markers/markers_NickName_OntologyName.csv"
      temp.table <- read.csv(fname,header=F);       temp.table <- as.matrix(temp.table)
      temp.matrix <- matrix(0,length(temp.table[,1])+1,2)
      temp.matrix[1:length(temp.table[,1]),] <- temp.table
      temp.table <- temp.matrix
      temp.table[length(temp.table[,1]),1] <- as.character(marker.list2[["Positive"]][q1]);    
      temp.table[length(temp.table[,1]),2] <-paste(temp.marker)
      write.table(temp.table, fname, sep = ",", col.names = FALSE, row.names = FALSE)  
        
      marker.list[["Positive"]][q1] <-temp.marker

    }
    if (nrow(res)>=2){
      cat("Marker", marker.list[["Positive"]][q1],"has not been changed, multiple possible markers in label\n")
    }
    if (nrow(res)==0 | nrow(res)>=2){
      temp.marker <- marker.list[["Positive"]][q1]
      # Define the cell.ctde.net SPARQL endpoint
      endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"
      # The options for 'read.csv' ensure that lines remain strings and not factors,
      # the lines are parsed based on the new line character (not a comma),
      # any quotes found in the parsing of the text are left intact,
      # and none of the lines being read should be treated as a header.
      que <- read.csv("Query/hasProperSynonym.txt", as.is = TRUE, sep = "\n", quote = "", 
                    header = FALSE)[, 1]
      prefix.info <- read.csv('Query/prefixInfo.txt', as.is = TRUE, sep = "\n",
                              quote = "", header = FALSE)[, 1]
      # Concatenate the query preceeded by all prefix information as a single string
      # to be passed to SPARQL
      query <- paste(c(prefix.info, que), collapse="\n")
      
      # Prepare marker for query by ensuring the marker is either followed by the  
      # end of the line or it has a symbol other than a letter or number after it, 
      # as that may indicate a different marker
      temp.marker <- paste(temp.marker, "$|", temp.marker, "[^0-9a-zA-Z-+][^ abes]", collapse="", sep="")
      query <- gsub("\\$marker", temp.marker, query)
      
      # Execute query
      res <- SPARQL(url=endpoint, query)$results
      
      # Small loop to check if the query is giving multiple label names. In this case the marker will not be changed
      if (nrow(res)>=1){
      temp = res[1,2]
      for(q5 in 1:nrow(res)){
        if(temp == res[q5,2]){
          temp = res[q5,2]
          sameLabels <- TRUE;
        }
        else{
          sameLabels <- FALSE;
          break
        }
      }
      
      if (sameLabels == TRUE){ # label exists and there is only one label
        temp.marker <- res[1,2]
        cat("Marker", marker.list[["Positive"]][q1], "has been updated to", temp.marker, "\n")
          
        # Update the markers_NickName_OntologyName.csv file
        fname <- "CL_Markers/markers_NickName_OntologyName.csv"
        temp.table <- read.csv(fname,header=F);       temp.table <- as.matrix(temp.table)
        temp.matrix <- matrix(0,length(temp.table[,1])+1,2)
        temp.matrix[1:length(temp.table[,1]),] <- temp.table
        temp.table <- temp.matrix
        temp.table[length(temp.table[,1]),1] <- as.character(marker.list2[["Positive"]][q1]);    
        temp.table[length(temp.table[,1]),2] <-paste(temp.marker)
        write.table(temp.table, fname, sep = ",", col.names = FALSE, row.names = FALSE)  
          
        marker.list[["Positive"]][q1] <-temp.marker
      }
      if (sameLabels == FALSE){ # label exists however there is two or more labels
        cat("Marker", marker.list[["Positive"]][q1],"has not been changed, multiple possible markers in synonyms\n")
      }
      }
      if (nrow(res)==0){ # label does not exist
          cat("Marker", marker.list[["Positive"]][q1],"could not be found\n")
      }
    }
  } # for loop
} # if statement

# The negative case
  if(length(marker.list[["Negative"]])!=0){
  for(q1 in 1:length(marker.list[["Negative"]])){
    
    if(file.exists(fname)){
    markers_NickName_OntologyName <- read.csv(fname , header=F)    
    markers_NickName_OntologyName <- as.matrix(markers_NickName_OntologyName)
 
    # skips query if there is a file named "markers_NickName_OntologyName" with all the marker ontology labels
    temp.location <- which(markers_NickName_OntologyName[,1]==marker.list2[["Negative"]][q1])
    if(length(temp.location)==1){
        cat("Marker", marker.list[["Negative"]][q1], "is called", markers_NickName_OntologyName[temp.location,2], "\n")
        marker.list[["Negative"]][q1] <- markers_NickName_OntologyName[temp.location,2]
        next
    }  
    }        
    # Change the short name markers in marker.list to the marker labels in the Ontology
    # This is only needed for the ones that the code has trouble finding
    # Hopefully with an updated Ontology these next three line can be deleted
    updatedException <- applyOntologyExceptionsNegative1(marker.list[["Negative"]], q1)
    marker.list[["Negative"]] <- applyOntologyExceptionsNegative2(marker.list[["Negative"]], q1)
    if ( updatedException == TRUE){next}
    
    temp.marker <- marker.list[["Negative"]][q1]
    # Define the cell.ctde.net SPARQL endpoint
    endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"
    
    # The options for 'read.csv' ensure that lines remain strings and not factors,
    # the lines are parsed based on the new line character (not a comma),
    # any quotes found in the parsing of the text are left intact,
    # and none of the lines being read should be treated as a header.
    que <- read.csv("Query/hasProperLabel.txt", as.is = TRUE, sep = "\n", quote = "", 
                  header = FALSE)[, 1]
    prefix.info <- read.csv('Query/prefixInfo.txt', as.is = TRUE, sep = "\n",
                            quote = "", header = FALSE)[, 1]
    # Concatenate the query preceeded by all prefix information as a single string
    # to be passed to SPARQL
    query <- paste(c(prefix.info, que), collapse="\n")
    
    # Prepare marker for query by ensuring the marker is either followed by the  
    # end of the line or it has a symbol other than a letter or number after it, 
    # as that may indicate a different marker
    temp.marker <- paste(temp.marker, "$|", temp.marker, "[^0-9a-zA-Z-+][^ abes]", collapse="", sep="")
    query <- gsub("\\$marker", temp.marker, query)
    
    # Execute query    
    res <- SPARQL(url=endpoint, query)$results
    
    if (nrow(res)==1){
       temp.marker <- res[1,2]
       cat("Marker", marker.list[["Negative"]][q1], "has been updated to", temp.marker, "\n")
    
       # Update the markers_NickName_OntologyName.csv file
       fname <- "CL_Markers/markers_NickName_OntologyName.csv"
       temp.table <- read.csv(fname,header=F);       temp.table <- as.matrix(temp.table)
       temp.matrix <- matrix(0,length(temp.table[,1])+1,2)
       temp.matrix[1:length(temp.table[,1]),] <- temp.table
       temp.table <- temp.matrix
       temp.table[length(temp.table[,1]),1] <- as.character(marker.list2[["Negative"]][q1]);    
       temp.table[length(temp.table[,1]),2] <- paste(temp.marker)
       write.table(temp.table, fname, sep = ",", col.names = FALSE, row.names = FALSE)          
        
       marker.list[["Negative"]][q1] <-temp.marker
    }
    if (nrow(res)>=2){
      cat("Marker", marker.list[["Negative"]][q1],"has not been changed, multiple possible markers in label\n")
    }
    
    if (nrow(res)==0|nrow(res)>=2){
      temp.marker <- marker.list[["Negative"]][q1]
      # Define the cell.ctde.net SPARQL endpoint
      endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"
      
      # The options for 'read.csv' ensure that lines remain strings and not factors,
      # the lines are parsed based on the new line character (not a comma),
      # any quotes found in the parsing of the text are left intact,
      # and none of the lines being read should be treated as a header.
      que <- read.csv("Query/hasProperSynonym.txt", as.is = TRUE, sep = "\n", quote = "", 
                    header = FALSE)[, 1]
      prefix.info <- read.csv('Query/prefixInfo.txt', as.is = TRUE, sep = "\n",
                              quote = "", header = FALSE)[, 1]
      # Concatenate the query preceeded by all prefix information as a single string
      # to be passed to SPARQL
      query <- paste(c(prefix.info, que), collapse="\n")
      
      # Prepare marker for query by ensuring the marker is either followed by the  
      # end of the line or it has a symbol other than a letter or number after it, 
      # as that may indicate a different marker
      temp.marker <- paste(temp.marker, "$|", temp.marker, "[^0-9a-zA-Z-+][^ abes]", collapse="", sep="")
      query <- gsub("\\$marker", temp.marker, query)
      
      # Execute query
      res <- SPARQL(url=endpoint, query)$results
      
      # small loop to check if the multiple hits are giving the same result
      if (nrow(res)>=1){
        temp = res[1,2]
        for(q5 in 1:nrow(res)){
          if(temp == res[q5,2]){
            temp = res[q5,2]
            sameLabels <- TRUE;
          }
          else{
            sameLabels <- FALSE;
            break
          }
        }
        if (sameLabels == TRUE){ # label exists and there is only one label
          temp.marker <- res[1,2]
          cat("Marker", marker.list[["Negative"]][q1], "has been updated to", temp.marker, "\n")
            
          # Update the markers_NickName_OntologyName.csv file
          fname <- "CL_Markers/markers_NickName_OntologyName.csv"
          temp.table <- read.csv(fname,header=F);       temp.table <- as.matrix(temp.table)
          temp.matrix <- matrix(0,length(temp.table[,1])+1,2)
          temp.matrix[1:length(temp.table[,1]),] <- temp.table
          temp.table <- temp.matrix
          temp.table[length(temp.table[,1]),1] <- as.character(marker.list2[["Negative"]][q1]);    
          temp.table[length(temp.table[,1]),2] <- paste(temp.marker)
          write.table(temp.table, fname, sep = ",", col.names = FALSE, row.names = FALSE)
            
          marker.list[["Negative"]][q1] <-temp.marker
        }
        if (sameLabels == FALSE){ # label exists however there is two or more labels
          cat("Marker", marker.list[["Negative"]][q1]," has not been changed, multiple possible markers in synonyms\n")
        }
      }
      if (nrow(res)==0){ # label does not exist
            cat("Marker", marker.list[["Negative"]][q1],"could not be found\n")                             
      }
    }
  }
  }

  return(marker.list)  
}

#################################################################################
# Similar to queryMarker, but for finding the parents for a specific cell population
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
  que <- read.csv(query.file, as.is = TRUE, sep = "\n", quote = "", 
                header = FALSE)[, 1]
  prefix.info <- read.csv('Query/prefixInfo.txt', as.is = TRUE, sep = "\n",
                          quote = "", header = FALSE)[, 1]
  # Concatenate the query preceeded by all prefix information as a single string
  # to be passed to SPARQL
  query <- paste(c(prefix.info, que), collapse="\n")
  # Add "^" to beginning and "$" to end to find an exact match for the label
  child.label <- paste("^", child.label, "$", sep = "", collapse = "")
  query <- gsub("\\$label", child.label, query)
  # Execute query
  res <- SPARQL(url=endpoint, query)$results 
  return (res)
}

#################################################################################
# Makes a query with SPARQL to the CL Ontology, and returns the results
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
  que <- read.csv(query.file, as.is = TRUE, sep = "\n", quote = "", 
                header = FALSE)[, 1]
  prefix.info <- read.csv('Query/prefixInfo.txt', as.is = TRUE, sep = "\n",
                          quote = "", header = FALSE)[, 1]
  # Concatenate the query preceeded by all prefix information as a single string
  # to be passed to SPARQL
  query <- paste(c(prefix.info, que), collapse="\n")

  # Prepare marker for query by ensuring the marker is either followed by the  
  # end of the line or it has a symbol other than a letter or number after it, 
  # as that may indicate a different marker
  marker <- paste(marker, "$|", marker, "[^0-9a-zA-Z-+][^ abes]", collapse="", sep="")
  query <- gsub("\\$marker", marker, query)
  if (!is.null(celllabel)){
    celllabel <- paste("^", celllabel, "$", sep = "", collapse = "")
    query <- gsub("\\$celllabel", celllabel, query)
  }

  # Execute query
  res <- SPARQL(url=endpoint, query)$results
  return (res)
}


