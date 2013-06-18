\name{Find.TSP.Path}
\alias{Find.TSP.Path}

\title{
Find.TSP.Path
}
\description{
Employs a heuristic approach to solve the traveling salesman problem.
}
\usage{
Find.TSP.Path(PositionList, mutation.matrix, insertion.type = "cheapest_insertion",
			  fix.start.pos = "Y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PositionList}{
		A dataframe consisting of six columns: 1) Residue Name, 2) Amino Acid number in the protein, 3) Side Chain, 4) X-coordinate, 5) Y-coordinate and 6) Z-coordinate. Please see \emph{get.Positions} and \emph{get.AlignedPositions} in the \emph{iPAC} package for further information on how to construct this matrix. 
}
  \item{mutation.matrix}{
		A matrix of 0's (no mutation) and 1's (mutation) where each column represents an amino acid in the protein and each row represents an individual sample (test subject, cell line, etc). Thus if column i in row j had a 1, that would mean that the ith amino acid for person j had a nonsynonomous mutation.
}
  \item{insertion.type}{
   Specifies the type of insertion method used. Please see the \emph{TSP} package for more details.
}
  \item{fix.start.pos}{
	The TSP package starts the path at a random amino acid. Such that the results are easily reproducible, the default starts the path on the first amino acid in the protein.
}
}
%\details{
% 
%}
\value{
\item{candidate.path}{A numeric vector of the sequence found through the protein.}
\item{candidate.path.distance}{The distance traveled along the candidate path.}
\item{dist.matrix}{The distance matrix between any two pairwise amino acids.}
\item{linear.path.distance}{The distance traveled if one were to visit the amino acids in the original sequence (1 -> 2 -> 3 -> ...->N}
}
\references{
Michael Hahsler and Kurt Hornik (2011). Traveling Salesperson Problem
(TSP) R package version 1.0-7.  \url{http://CRAN.R-project.org/}.

Gregory Ryslik and Hongyu Zhao (2012). iPAC: Identification of Protein Amino acid Clustering. R package version 1.1.3.  \url{http://www.bioconductor.org/}.
}

\examples{
#Load the position and mutational data
CIF<-"http://www.pdb.org/pdb/files/3GFT.cif"
Fasta<-"http://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.Positions<-get.Positions(CIF,Fasta, "A")
data(KRAS.Mutations)

#Save all the results to path.results
path.results <- Find.TSP.Path(KRAS.Positions$Positions, KRAS.Mutations)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Graph }
\keyword{ Traveling Salesman }
