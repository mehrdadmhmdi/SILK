\name{theme_silk}
\alias{theme_silk}
\alias{method_palette}
\alias{save_plot}
\alias{make_prediction_plots}
\alias{make_dgm_plots}
\title{SILK Visualization Functions}
\usage{
theme_silk(base_size = NULL)
method_palette()
save_plot(plot, out_dir, name, width = 12, height = 7)
make_prediction_plots(out_dir)
make_dgm_plots(out_dir, src_dir)
}
\arguments{
  \item{base_size}{Numeric. Base font size for theme.}
  \item{plot}{A ggplot object.}
  \item{out_dir}{Character. Output directory.}
  \item{name}{Character. File name without extension.}
  \item{width}{Numeric. Figure width in inches.}
  \item{height}{Numeric. Figure height in inches.}
  \item{src_dir}{Character. Directory for source data CSVs.}
}
\value{Various; see individual function descriptions.}
\description{Publication-quality visualization functions for SILK results. Requires ggplot2 (in Suggests).}
