% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_free_energy_change}
\alias{plot_free_energy_change}
\title{Plot change in free energy over time}
\usage{
plot_free_energy_change(results)
}
\arguments{
\item{results}{list; object containing result  (output of read_output)}
}
\value{
A ggplot2 object
}
\description{
\code{plot_free_energy_change} Plot change in free energy over time
}
\examples{
data(results)
plot_free_energy_change(results)

}
