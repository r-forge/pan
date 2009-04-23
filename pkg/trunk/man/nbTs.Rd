\name{nbTs}
\alias{nbTs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ nbTs }
\description{
  Returns the number or the proportion of transitions in a DNA sequence alignment.
}
\usage{
nbTs(x, pairwise.deletion = FALSE, as.matrix = FALSE, scaled = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a matrix or a list containing DNA sequences. }
  \item{pairwise.deletion}{ a logical indicating whether to delete the sites
          with missing data in a pairwise way. The default is to delete
          the sites with at least one missing data for all sequences. }
  \item{as.matrix}{ a logical indicating whether to return the results as a
          matrix. The default is to return an object of class dist.}
  \item{scaled}{ a logical indicating whether to return the proportion of sites with
                 transitions or the number of transitions. }
}

\value{
  An object of class \code{dist} or a numeric matrix if \code{as.matrix = TRUE}.
}

\author{ Francois Michonneau <francois.michonneau@gmail.com> modified from
         code written by Emmanuel Paradis. }

\seealso{ \code{\link{dist.dna}} }
\examples{
        xx <- 10

}
