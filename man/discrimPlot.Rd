\name{discrimPlot}
\alias{discrimPlot}
%\alias{d.primePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensory discrimination power analysis}
\description{
  Plots the power and the sample size of a difference or similarity test for a sensory
  discrimination experiment using the binomial distribution.

}
\usage{
discrimPlot(type=c("power-SS","power-ES","SS-ES"),display=c("protocol","alpha","ES","power","SS"),
            pdA=c(0.25,0.375,0.5), pd0 = 0, target.power = c(0.9,0.85,0.8), sample.size=c(20,50,100),  
            alpha = c(0.05,0.10,0.01), pGuess = c(1/3,1/2,1/9,1/4), test = c("difference", "similarity"), 
            statistic = c("exact","stable.exact", "normal", "cont.normal"))


}
\arguments{
\item{type}{the type of plot: "power-SS" corresponds to power vs. sample size, "power-ES" to power vs effect size (pdA) and "SS-ES" to sample size vs. effect size (pdA)}
\item{display}{the range to display on the graph: "protocol" displays the probability of guessing of different protocols, "alpha" different type I levels, "ES" different levels of effect size (pdA), "power" different levels of test powers and "SS" different sample size}
  \item{pdA}{the probability of discrimination for the model under the
    alternative hypothesis; scalar between zero and one}
  \item{pd0}{the probability of discrimination under the
    null hypothesis; scalar between zero and one}
  \item{target.power}{the desired power for the test}
  \item{sample.size}{the sample size; a scalar positive integer}
  \item{alpha}{the type I level of the test; scalar between zero and
    one}
  \item{pGuess}{the guessing probability for the discrimination
    protocol, e.g. 1/2 for duo-trio and 2-AFC, 1/3 for
    triangle, tetrad and 3-AFC, 1/4 for the double duo-trio and the double 2-AFC and 1/9 for the double
    triangle, double tetrad and double 3-AFC; scalar between zero and one}
  \item{test}{the type of one-sided binomial test (direction of the
    alternative hypothesis): "difference" corresponds "greater" and
    "similarity" corresponds to "less"}
  \item{statistic}{should power determination be based on the 'exact'
           binomial test, the normal approximation to this, or the
           normal approximation with continuity correction?}
}
\details{
  The power of the standard one-tailed difference test where the null
  hypothesis is "no difference" is obtained with \code{pd0 = 0}.

  For difference tests \code{pdA} or \code{d.primeA} (the sensory
  difference under the alternative hypothesis) has to be larger than
  \code{pd0} or \code{d.prime0} (the sensory  difference under the null
  hypothesis). The sample size of the standard one-tailed difference
  test where the null hypothesis of "no difference" is obtained with
  \code{pd0 = 0} or \code{d.prime0 = 0}.

  For similarity tests it is required that \code{pd0} > \code{pdA} or
  equivalently that \code{d.prime0} > \code{d.primeA}. Here, the
  interval [0, \code{pdA}] or [0, \code{d.primeA}]
  is the similarity region covering sensory differences for which we
  would say that the products are similar.

  The probability of a correct answer under the null hypothesis is
  given by \code{pd0 + pGuess * (1 - pd0)}. Similarly, the probability
  of a correct answer under the alternative hypothesis is given by
  \code{pdA + pGuess * (1 - pdA)}.

  The \code{statistic} argument:
  \itemize{
    \item{\code{"exact"}}{
      is the conventional sample size for the exact binomial test:
      The smallest sample size that gives the desired power
      (\code{target.power}) at the given
      significance level. Ususally slightly higher sample sizes will not
      have the desired power, however. This is due to the non-monotonic
      behavior of power as a function of sample size.
    }
    \item{\code{"stable.exact"}}{
      is so-called stable exact sample size proposed by Ennis and
      Jesionka (2011) which has the property that no larger sample sizes
      has a power less than the \code{target.power}.
    }
    \item{\code{"normal"}}{
      is the normal approximation to the exact binomial sample size
      without any continuity adjustment. This usually provides a sample
      size that is smaller than the sample size for the exact binomial
      test.
    }
    \item{\code{"cont.normal"}}{
      is the continuity adjusted normal approximation to
      the sample size for the exact binomial test. This sample size is
      usually closer to the exact sample size than the unadjusted
      approximation and usually higher than the unadjusted
      approximation.
    }
  }
}
\value{
  An interactive graph.
}
\references{
  Brockhoff, P.B. and Christensen, R.H.B (2010). Thurstonian
  models for sensory discrimination tests as generalized linear models.
  Food Quality and Preference, 21, pp. 330-338.
  
  Ennis, J.M. and V. Jesionka (2011). The power of sensory
  discrimination methods revisited. Journal of Sensory Studies, 26,
  pp. 371-382.
}
\author{Sophie Birot}

\seealso{\code{\link{findcr}}, \code{\link{discrimPwr}},
  \code{\link{discrim}}, \code{\link{discrimSim}},
  \code{\link{samediff}},
  \code{\link{AnotA}}, \code{\link{discrimSS}}
}

\examples{
## Plotting power vs. sample size for different pGuess
discrimPlot(type="power-SS",display="protocol")

## Plotting sample size vs. effect size for different alpha with stable.exact statistic
discrimPlot(type="SS-ES",display="alpha",statistic = "stable.exact")

}
\keyword{models}% at least one, from doc/KEYWORDS
