\name{GraphClust}
\alias{GraphClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GraphClust
}
\description{
Finds mutational clusters after reordering the protein using the traveling salesman approach.
}
\usage{
GraphClust(mutation.data, position.data, insertion.type = "cheapest_insertion", alpha = 0.05,
		   MultComp = "Bonferroni", fix.start.pos = "Y", Include.Culled = "Y",
		   Include.Full = "Y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutation.data}{
 A matrix of 0's (no mutation) and 1's (mutation) where each column represents an amino acid in the protein and each row represents an individual sample (test subject, cell line, etc). Thus if column i in row j had a 1, that would mean that the ith amino acid for person j had a nonsynonomous mutation.
}
  \item{position.data}{
A dataframe consisting of six columns: 1) Residue Name, 2) Amino Acid number in the protein, 3) Side Chain, 4) X-coordinate, 5) Y-coordinate and 6) Z-coordinate. Please see \emph{get.Positions} and \emph{get.AlignedPositions} in the \emph{iPAC} package for further information on how to construct this matrix.
}
  \item{insertion.type}{
 Specifies the type of insertion method used. Please see the \emph{TSP} package for more details.
}
  \item{alpha}{
The significance level required in order to find a mutational cluster significance. Please see the NMC package for further information.
}
  \item{MultComp}{
The multiple comparison adjustment required as all pairwise mutations are considered. Options are: ``Bonferroni", "BH", or "None".
}
  \item{fix.start.pos}{
The TSP package starts the path at a random amino acid. Such that the results are easily reproducible, the default starts the path on the first amino acid in the protein.
}
  \item{Include.Culled}{
If "Y", the standard NMC algorithm will be run on the protein after removing the amino acids for which there is no positional data.
}
  \item{Include.Full}{
If "Y", the standard NMC algorithm will be run on the full protein sequence.
}
}
\details{
The protein reordering is done using the \emph{TSP} package available on CRAN. This hamiltonian path then serves as the new protein ordering.

The position data can be created via the ``get.AlignedPositions" or the ``get.Positions" functions available via the imported \emph{iPAC} package.

The mutation matrix must have the default R column headings ``V1", ``V2",...,``VN", where N is the last amino acid in the protein. No positions should be skipped in the mutaion matrix.

When unmapping back to the original space, the end points of the cluster in the mapped space are used as the endpoints of the cluster in the unmapped space.
}

\value{
\item{Remapped}{This shows the clusters found while taking the 3D structure into account and remapping the protein using a traveling salesman approach.}
\item{OriginalCulled}{This shows the clusters found if you run the NMC algorithm on the canonical linear protein, but with the amino acids for which we don't have 3D positional data removed.}
\item{Original}{This shows the clusters found if you run the NMC algorithn on the canonical linear protein with all the amino acids.}
\item{candidate.path}{This shows the path found by the TSP package that heuristically minimizes the total distance through the protein.}
\item{path.distance}{The length of the candidate path if traveled from start to finish.}
\item{linear.path.distance}{The length of the sequential path 1,2,3...,N (where N is the total number of amino acids in the protein).}
\item{protein.graph}{A graph object created by the \emph{igraph} package that has edges between amino acids on the candidate.path. This can be passed to plotting functions to create visual represnetations.}
\item{missing.positions}{This shows which amino acids are present in the mutation matrix but for which we do not have positions. These amino acids are cut from the protein when calculating the \emph{Remapped} and \emph{OriginalCulled} results.}
}
\references{
Ye et. al., Statistical method on nonrandom clustering with application to somatic mutations in cancer. \emph{BMC Bioinformatics}. 2010. doi:10.1186/1471-2105-11-11.

Michael Hahsler and Kurt Hornik (2011). Traveling Salesperson Problem (TSP) R package version 1.0-7.  \url{http://CRAN.R-project.org/}.

Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. \url{http://igraph.sf.net} 

Gregory Ryslik and Hongyu Zhao (2012). iPAC: Identification of Protein Amino acid Clustering. R package version 1.1.3. \url{http://www.bioconductor.org/}.
}

\examples{
\dontrun{
#Load the positional and mutatioanl data
CIF<-"http://www.pdb.org/pdb/files/3GFT.cif"
Fasta<-"http://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.Positions<-get.Positions(CIF,Fasta, "A")
data(KRAS.Mutations)

#Calculate the required clusters
GraphClust(KRAS.Mutations,KRAS.Positions$Positions,insertion.type = "cheapest_insertion",
	   alpha = 0.05, MultComp = "Bonferroni")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clusters }
\keyword{ Mutations }% __ONLY ONE__ keyword per line
