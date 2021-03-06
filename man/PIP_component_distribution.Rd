% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{PIP_component_distribution}
\alias{PIP_component_distribution}
\title{Plot PIP distribution for a single component}
\usage{
PIP_component_distribution(results, component, omic = 1)
}
\arguments{
\item{results}{list; object containing result  (output of read_output)}

\item{component}{integer or character; index or name of the component you want to plot}

\item{omic}{integer or character; index of the omic type to calculate contributions for, or name of omic type if you assigned a name}
}
\value{
A ggplot2 object
}
\description{
\code{PIP_component_distribution} Plot final PIP distribution for a single component
}
\examples{
data(results)
PIP_component_distribution(results, 8)

}
