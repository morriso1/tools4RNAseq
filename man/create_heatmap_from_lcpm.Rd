% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAseq_plotting.R
\name{create_heatmap_from_lcpm}
\alias{create_heatmap_from_lcpm}
\title{Create a heatmap using the heatmap.2 function from gplots}
\usage{
create_heatmap_from_lcpm(
  loc_exp_lcpm,
  groups,
  sam_ids,
  gene_set_charV,
  reg_express = "LacZRNAi",
  col_cell = NA,
  ...
)
}
\arguments{
\item{loc_exp_lcpm}{Log transformed count matrix.}

\item{groups}{Column of DGElist$samples that contains group info.}

\item{sam_ids}{Column of DGElist$samples that contains sam_id info.}

\item{gene_set_charV}{Gene set as a vector. Identifies need to match rownames of loc_exp_lcpm e.g. flybase_ids}

\item{reg_express}{REGEX expression to identify columns of groups to include in heatmap}

\item{col_cell}{Vector of colors to label each of the sam_id columns in the heatmap}

\item{...}{Passing the dots to heatmap.2 function.}
}
\value{
plot
}
\description{
Create a heatmap using the heatmap.2 function from gplots
}
