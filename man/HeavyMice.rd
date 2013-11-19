\name{HeavyMice}
\alias{HeavyMice}
\docType{data}
\title{Heavy Mice Dataset.}
\description{This data set gives a subset of the microarray expression data from
the liver tissue of heavy female mice.  There were 3421 genes and 135 mice in the full data set; this data set was obtained by removing genes and mice with missing values.  Then the 50 heaviest mice (weights greater than 40.5 g) were selected.  Finally, univariate regressions of mouse weights on each individual gene were performed and 314 genes with z-scores greater than 5 were selected.}
\usage{data(HeavyMice)}
\format{A 50 by 314 matrix.}
\source{
The complete mouse weight data set can be found at 
\url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/MouseWeight/}.
}
\references{  
Fuller, T., Ghazalpour, A., Aten, J., Drake, T., Lusis, A., and Horvath, S. (2007) Weighted gene coexpression network analysis strategies applied to mouse weight. \emph{Mammalian Genome}, \bold{28}, 463--472. 

Ghazalpour, A., Doss, S., Zhang, B., Wang, S., Plaisier, C., Castellanos, R., Brozell, A., Schadt, E.E., Drake, T.A., Lusis, A.J., and Horvath, S. (2006) Integrated genetic and network analysis to characterize genes related to mouse weight. \emph{PLoS Genetics}, \bold{2}, 130.

Gill, R., Datta, S., and Datta, S. (2010) A statistical framework for differential network analysis from microarray data. \emph{BMC Bioinformatics}, \bold{11}, 95.
}
\keyword{datasets}
