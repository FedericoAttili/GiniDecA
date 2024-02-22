\name{GiniDecA-package}
\alias{GiniDecA-package}
\alias{GiniDecA}
\docType{package}
\title{
  A short title line describing what the package does
}
\description{
  A more detailed description of what the package does. A length
  of about one to five lines is recommended.
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
Federico Attili
Maintainer: Federico Attili <federico.attili2@unibo.it>
}
\references{
  This optional section can contain literature or other references for
  background information.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\usage{
\dontrun{giniDec(x, groups, contrib = T)}}
\examples{
\dontrun{
## Example usage

# Generating 4 random numbers between 1 and 10
x = runif(4, 1, 10)

# Creating a vector of group identifiers (1, 2, 1, 2)
groups = c(rep(c(1:2), 2))


# Computing the Gini decomposition
res1 = giniDec(x, groups, contrib = T)


##   Example with different-sized groups and weights

# Generating 5 random numbers between 1 and 10
x = runif(5, 1, 10)

# Creating a vector of group identifiers (1, 2, 1, 2, 2)
groups = c(rep(c(1:2), 2), 2)

# Creating a weight vector (1, 2, 1, 2, 2)
w = c(rep.int(c(1, 2), 2), 2)

# Computing the Gini decomposition with additional parameters
res2 = giniDec(x, groups, w = w, n_equalizer = 'quant', type = 7, contrib = T)

}
}
