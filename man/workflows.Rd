\name{sparq_assess_workflow}
\alias{sparq_assess_workflow}
\alias{sparq_preset}
\alias{sparq_validate_workflow}
\alias{sparq_assess_adaptive}
\alias{sparq_assess_script}
\alias{sparq_prepare_spatial}
\alias{sparq_subset_spatial}
\alias{sparq_assess_spatial}
\title{Run SPARQ through workflow, script, adaptive, and spatial-object interfaces}
\description{
These high-level interfaces run a reference analysis, perturb spatial observations,
rerun the supplied analysis, compare the outputs, and return a standard SPARQ
assessment. They do not alter the original input object.
}
\section{Calling conventions}{\preformatted{
sparq_assess_workflow(data, analysis_function,
  output_type = c("scalar", "ranked", "partition", "graph"),
  comparator = NULL, ...)

sparq_preset(preset = c("quick", "standard", "thorough"))

sparq_validate_workflow(data, analysis_function,
  output_type = c("scalar", "ranked", "partition", "graph"), ...)

sparq_assess_adaptive(data, analysis_function, comparator,
  support_threshold = NULL, reliability_tolerance = NULL, ...)

sparq_assess_script(data, script, input_object, result_object,
  output_type = c("scalar", "ranked", "partition", "graph"), ...)

sparq_prepare_spatial(object, object_type = "auto", coordinate_columns = NULL,
  spot_id_column = NULL, image = NULL, adapter = NULL)

sparq_subset_spatial(spatial_input, spot_ids)

sparq_assess_spatial(object, analysis_function,
  output_type = c("scalar", "ranked", "partition", "graph"), ...)
}}
\arguments{
  \item{data}{A nonempty spot-level data frame. Each row is an observation that
  can be retained or removed.}
  \item{analysis_function}{A function accepting the complete or perturbed input
  and returning the selected output type.}
  \item{output_type}{The result representation: one numeric scalar, an ordered
  feature vector, named spot labels, or canonical graph edge IDs.}
  \item{preset}{One of \code{"quick"} (25 perturbations), \code{"standard"}
  (100 perturbations), or \code{"thorough"} (500 perturbations). A preset sets
  retention, minimum iterations, bootstrap budgets, and progress frequency;
  explicit arguments override an individual preset value.}
  \item{preflight_retention}{The fraction retained in the one small dry-run
  perturbation; by default 95 percent.}
  \item{comparator}{A function comparing reference and perturbed results.}
  \item{support_threshold}{A prespecified decision threshold for adaptive
  stopping. It is never inferred from the reference result.}
  \item{reliability_tolerance}{Optional prespecified normalized comparator
  distance defining an acceptably preserved result. When supplied to an
  assessment, SPARQ returns a reliability decision. When supplied to adaptive
  assessment, it replaces support-score stopping with anytime-valid
  beta-binomial-mixture reliability stopping.}
  \item{required_reliability}{Minimum preservation probability required for a
  tolerance-based reliability decision.}
  \item{reliability_confidence}{Confidence level for the reliability interval.}
  \item{loss_model}{Optional explicit \code{sparq_loss_model} created with
  \code{sparq_define_loss_model}. It determines the retained-spot sampling
  distribution and is recorded with the assessment.}
  \item{script}{One existing R script.}
  \item{input_object}{The name under which SPARQ injects current perturbed data
  into a script environment.}
  \item{result_object}{The name of the result created by the script.}
  \item{object}{A data frame, Seurat, SpatialExperiment, Giotto, or custom
  spatial object.}
  \item{object_type}{Object type, or \code{"auto"} to detect it.}
  \item{coordinate_columns}{Optional names of x and y coordinate columns.}
  \item{spot_id_column}{Optional spot-ID column for a data frame.}
  \item{image}{Spatial image name for an object with multiple Seurat images.}
  \item{adapter}{For another object class, a list with \code{extract(object)}
  returning \code{spot_ids} and row-named \code{coordinates}, plus
  \code{subset(object, spot_ids)} returning the retained object.}
  \item{spatial_input}{An object returned by \code{sparq_prepare_spatial}.}
  \item{spot_ids}{Unique retained spot IDs.}
  \item{...}{Assessment controls including \code{retention},
  \code{n_iterations}, \code{stress_model}, \code{seed}, and \code{verbose}.}
}
\value{
\code{sparq_assess_workflow}, \code{sparq_assess_script}, and
\code{sparq_assess_spatial} return a SPARQ assessment with a reference result,
comparison table, support summary, failed-iteration ledger, and settings.
When a reliability tolerance is supplied, they also return a three-way
tolerance-based reliability decision without altering score calibration.
\code{sparq_assess_adaptive} also returns the stopping history, decision, and,
for reliability stopping, the full anytime-valid confidence sequence.
\code{sparq_preset} returns the exact named settings selected by a preset.
\code{sparq_validate_workflow} returns a \code{sparq_workflow_preflight}
object containing the reference and dry-run results plus a pass/fail table.
\code{sparq_prepare_spatial} returns a \code{sparq_spatial_input} that retains
the original object and a standardized spot-coordinate table.
}
\details{
The spatial adapter subsets the original object on each perturbation before
calling \code{analysis_function}; it does not convert a Seurat,
SpatialExperiment, or Giotto object into a data frame for the analysis itself.
Seurat requires \pkg{SeuratObject}; SpatialExperiment requires
\pkg{SpatialExperiment}; Giotto requires \pkg{GiottoClass} or \pkg{Giotto}.
These are optional dependencies, required only for their respective object type.
}
\examples{
spots <- data.frame(value = c(2, 3, 5, 7, 11, 13), x = 1:6, y = 6:1)

assessment <- sparq_assess_workflow(
  data = spots,
  analysis_function = function(x) mean(x$value),
  output_type = "scalar",
  retention = 0.75,
  n_iterations = 4,
  min_iterations = 3,
  instability_bootstrap_B = 20,
  verbose = FALSE
)
assessment

sparq_preset("quick")

preflight <- sparq_validate_workflow(
  data = spots,
  analysis_function = function(x) mean(x$value),
  output_type = "scalar",
  verbose = FALSE
)
preflight

spatial_input <- sparq_prepare_spatial(spots, coordinate_columns = c("x", "y"))
subsetted_spots <- sparq_subset_spatial(spatial_input, c("1", "2", "3"))
}
