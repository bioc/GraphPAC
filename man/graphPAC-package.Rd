\name{GraphPAC-package}
\alias{GraphPAC-package}
\alias{GraphPAC}
\docType{package}
\title{
Using Graph Theory to Identify Mutational Clusters of Amino Acids.
}
\description{
The \emph{GraphPAC} package identifies statistically significant clusters of non-synonomous amino acid mutations and is a sister package to \emph{iPAC}. \emph{GraphPAC} reorders the protein into a one dimensional space via a graph theoretrical approach. Specifically, the traveling salesman problem (TSP) is solved heuristically via the \emph{TSP} package. Once solved, the mutational data is reordered to follow the hamiltonian path and the nmc algorithm is run to find the mutational clusters on the remapped protein. Unlike the MDS remapping approach that is used in \emph{iPAC}, distant amino acids no longer have an effect on each other's position in one dimensional space allowing for a closer representation of the underlying biology.
}
\details{
Please see the documentation for ``get.Positions", ``get.AlignedPositions", and "Plot.Protein.Linear" in the iPAC package. There you will find information on getting basic positional data and plotting functions.
}
\author{
Gregory Ryslik Hongyu Zhao

Maintainer: Gregory A. Ryslik <gregory.ryslik@yale.edu>
}
\references{

Ye et. al., Statistical method on nonrandom clustering with application to somatic mutations in cancer. \emph{BMC Bioinformatics}. 2010. doi:10.1186/1471-2105-11-11.

Michael Hahsler and Kurt Hornik (2011). Traveling Salesperson Problem (TSP) R package version 1.0-7.  \url{http://CRAN.R-project.org/}.

Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. \url{http://igraph.sf.net}

Gregory Ryslik and Hongyu Zhao (2012). iPAC: Identification of Protein
Amino acid Clustering. R package version 1.1.3. \url{http://www.bioconductor.org/}.

Bioconductor: Open software development for computational biology and bioinformatics R. Gentleman, V. J. Carey, D. M. Bates, B.Bolstad, M. Dettling, S.
  Dudoit, B. Ellis, L. Gautier, Y. Ge, and others 2004, Genome Biology, Vol. 5, R80


}

\examples{
\dontrun{
#Load the positional and mutatioanl data
CIF<-"https://files.rcsb.org/view/3GFT.cif"
Fasta<-"https://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.Positions<-get.Positions(CIF,Fasta, "A")
data(KRAS.Mutations)

#Calculate the required clusters
GraphClust(KRAS.Mutations,KRAS.Positions$Positions,insertion.type = "cheapest_insertion",
		   alpha = 0.05, MultComp = "Bonferroni")
}
}
