\name{Plot.Protein}
\alias{Plot.Protein}

\title{
Plot.Protein
}
\description{
Creates a circular interactive plot of the path through the protein.
}
\usage{
Plot.Protein(graph, path, vertex.size = 5, color.palette = "heat")
}

\arguments{
  \item{graph}{
The graph object returned by GraphClust ($protein.graph).
}
  \item{path}{
The path returned by GraphClust ($candidate.path).
}
  \item{vertex.size}{
How large you want each vertex to be.
}
  \item{color.palette}{
Possible options are: "heat", "gray", "topo", "cm".
}
}
\details{
This will plot the amino acids in a circular directed graph. The vertices can be dragged around to enhance the visual representation. This is meant to complement the \emph{Plot.Protein.Linear} function in \emph{iPAC} which is also applicable in this package. 
}

\note{
This function is based on the ``tkplot" function in \emph{igraph}. Please see the documentation for that package for the necessary requirements. Special thanks to Dr. G\'{a}bor Cs\'{a}rdi (creator of the \emph{igraph} package) for his help.
}

\references{
Gregory Ryslik and Hongyu Zhao (2012). iPAC: Identification of Protein
Amino acid Clustering. R package version
1.1.3. \url{http://www.bioconductor.org/}.

Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. \url{http://igraph.sf.net}.
}

\examples{
\dontrun{
#Loads the mutational and positional data
CIF<-"https://files.rcsb.org/view/3GFT.cif"
Fasta<-"http://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.Positions<-get.Positions(CIF,Fasta, "A")
data(KRAS.Mutations)

#gets the cluster results and graph object
my.graph.clusters <- GraphClust(KRAS.Mutations,KRAS.Positions$Positions,
								insertion.type = "cheapest_insertion",alpha = 0.05,
								MultComp = "Bonferroni")

Plot.Protein(my.graph.clusters$protein.graph, my.graph.clusters$candidate.path,
			 vertex.size=5, color.palette="heat")
}
}

\keyword{ Order }
\keyword{ Graph }
\keyword{ Amino Acids }
