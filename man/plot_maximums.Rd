% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_maximums}
\alias{plot_maximums}
\title{Plot maximum score and loading by component}
\usage{
plot_maximums(results, omic = 1, labels = T, ...)
}
\arguments{
\item{results}{list; object containing result  (output of read_output)}

\item{omic}{numeric or character; index of the omic type to plot loading for, or name of omic type if you assigned a name}
}
\value{
A ggplot2 object
}
\description{
\code{plot_maximums} Plot maximum score and loading by component
}
\examples{
data(results)
plot_maximums(results)

}
