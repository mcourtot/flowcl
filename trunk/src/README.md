# Description of code scripts.
# Last updated: July 19, 2013 by Radina Droumeva (radina.droumeva@gmail.com)

Note: as a first step you must have R installed, as well as the R package "SPARQL":
http://cran.r-project.org/web/packages/SPARQL/

Once you install R on your system, run it (type R in the terminal in linux or
double click the icon in Windows/Mac) and try typing install.packages('SPARQL') 
into the terminal/command input and see if your system will automatically 
install it for you. 

If it does not, you will have to download the appropriate file from the link above
 and install it manually. You may need to use google to get instructions on how 
to install R packages on your operating system. Also note that this package depends 
on the R packages "XML", "bitops" and "RCurl". If the automatic installation 
(install.packages("SPARQL")) method above did not work, depending on the R version 
you have you may have to repeat the manual installation of the above packages in 
a sequential manner (bitops, RCurl, XML, SPARQL). You can test that the installation
was successful by opening R and typing "library('SPARQL')" and ensure no error
occurs (warnings are not errors).


Main script:

********************************************************************************
markerQueryTest.R: R script to test functionality. Near the top of the script,
you can change the variable "phenotye". This is the script which could eventually
evolve into an R package.
********************************************************************************


Collection of supporting functions:

********************************************************************************
rQueryFunctions.R: A collection of R-based SPARQL query functions:
  - queryMarker(marker, query.file, celllabel) fetches the unique (ID, label) 
    pairs with exact synonym matching the passed 'marker' (e.g. marker can be "CD19"),
    or checking for a single population with label 'celllabel' whether it lacks or
    has plasma membrane part some of marker 'm'.
  - parentQuery(child.label, query.file): returns the parents of a cell population
    'child.label'.

supportingFunctions.R: A collection of non-SPARQL R supporting functions:
  - phenoParse(phenotype) takes in a phenotype string such as "CD19+CD20-" and
    returns a list with two vectors -- one vector contains the markers with
    'positive' expression, and the other contains the markers with 'negative'
    expression. These markers can later be passed to 'queryMarker(marker)' above.
  - tabulateResults(res) organizes the initial queries to get populations with
    markers matching our phenotype of interest.
  - col2hex(colour, alpha) is a small helper function to make transparent colours,
    used for visualization purposes.
  - getPeaks(x, span) helper function in identifying a threshold on the population label
    scores. It locates the peaks in the density of scores.
  - getScoreCutoff(scores) uses getPeaks to find the cutoff for scores. The highest
    two peaks are located and then the minimum point in the density between them is
    used as the cutoff.
  - treePlot(...): this function creates the tree-structure plot that is the final
    result. This step should be improved/completely replaced with a different tool
    to visualize the final results.
********************************************************************************


Collection of SPARQL queries, all saved individually in .txt files to allow for
independent editting:

********************************************************************************
prefixInfo.txt: this is pasted into every query before the query is sent. Contains
    all prefix declarations used in the queries.
 
getMatchingSynonyms.txt: Contains a SPARQL query which is parsed into 
    'queryMarker(marker, ...)'. Kept separate in case it needs to be edited for 
    namespace issues or graph issues, etc.

getParentClasses.txt: Filters on cell label and returns its parents.

hasPlasmaMembranePart.txt: Filters on a marker label and finds populations which
    have plasma membrane part that marker.
    
lacksPlasmaMembranePart.txt: similarly to above...

hasPMPsingle.txt: Filters on marker label AND cell label, basically checks if
    the population of interest has plasma membrane part the marker of interest.
    
lacksPMPsingle.txt: similarly, but checks if it lacks the marker (explicitly lacks).
********************************************************************************


