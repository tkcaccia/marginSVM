#' List real spatial benchmark datasets
#'
#' Reports which real datasets can be redistributed with the package and records
#' the status of the other biological studies used by the publication scripts.
#'
#' @return A data frame with source, license, inclusion, and benchmark details.
#' @export
#' @examples
#' available_spatial_benchmarks()
available_spatial_benchmarks <- function() {
  data.frame(
    dataset = c("dlpfc", "merfish", "crc"),
    included = c(TRUE, FALSE, TRUE),
    observations = c(47329L, 28317L, 194541L),
    classes = c(7L, 8L, 19L),
    scenarios = c(45L, 45L, 60L),
    license = c(
      "Artistic-2.0 (spatialLIBD data package)",
      "CC0 raw Dryad data; no explicit license for derived BASS domain labels",
      "CC BY 4.0 (10x Genomics source and author-derived annotations)"
    ),
    source = c(
      "https://bioconductor.org/packages/spatialLIBD",
      "https://doi.org/10.5061/dryad.8t8s248",
      paste0(
        "https://www.10xgenomics.com/datasets/",
        "visium-hd-cytassist-gene-expression-libraries-of-human-crc-v4"
      )
    ),
    note = c(
      "Coordinates, layer labels, and frozen corruption inputs are bundled.",
      "Use the publication script with a local processed annotation file.",
      paste(
        "Coordinates, 19 WSI annotation labels, and deterministic corruption",
        "recipes are bundled; counts and tissue imagery are excluded."
      )
    ),
    stringsAsFactors = FALSE
  )
}

#' Load a bundled real spatial benchmark
#'
#' Loads one frozen corruption scenario from either the human dorsolateral
#' prefrontal cortex (DLPFC) Visium benchmark or the colorectal cancer (CRC)
#' Visium HD benchmark. The package contains coordinates and labels, not
#' expression counts or histology images.
#'
#' @param name Dataset name. `"dlpfc"` and `"crc"` are bundled. The legacy
#'   `"colorectal"` alias is accepted for compatibility.
#' @param scenario Scenario number or scenario identifier. Inspect the `design`
#'   element of the corresponding data object for the complete design.
#'
#' @return A `spatial_refinement_benchmark` ready for
#'   [benchmark_spatial_refiners()].
#' @export
#' @examples
#' dlpfc <- load_spatial_benchmark("dlpfc", scenario = 1)
#' dim(dlpfc$xy)
#' \donttest{
#' crc <- load_spatial_benchmark("crc", "CRC_random_25_r1")
#' mean(crc$labels != crc$truth)
#' }
load_spatial_benchmark <- function(
    name = c("dlpfc", "merfish", "crc", "colorectal"), scenario = 1L) {
  requested_name <- match.arg(name)
  name <- if (identical(requested_name, "colorectal")) "crc" else requested_name
  if (name == "merfish") {
    status <- available_spatial_benchmarks()
    row <- status[status$dataset == name, ]
    stop("`", name, "` is not bundled: ", row$license, ". ", row$note,
         call. = FALSE)
  }

  environment <- new.env(parent = emptyenv())
  object_key <- if (name == "crc") "colorectal" else name
  object_name <- paste0(object_key, "_benchmark")
  utils::data(list = object_name, package = "fibermargin", envir = environment)
  if (!exists(object_name, envir = environment, inherits = FALSE)) {
    stop("The bundled ", name, " benchmark could not be loaded.", call. = FALSE)
  }
  source <- get(object_name, envir = environment, inherits = FALSE)
  if (length(scenario) != 1L || is.na(scenario)) {
    stop("`scenario` must be one scenario number or identifier.", call. = FALSE)
  }
  scenario_index <- if (is.character(scenario)) {
    match(scenario, source$design$scenario_id)
  } else {
    value <- as.integer(scenario)
    if (!is.na(value) && value >= 1L && value <= nrow(source$design)) value else NA_integer_
  }
  if (is.na(scenario_index)) {
    stop("Unknown ", name, " benchmark scenario. Use a number from 1 to ",
         nrow(source$design), " or a value from `design$scenario_id`.",
         call. = FALSE)
  }

  scenario_design <- source$design[scenario_index, , drop = FALSE]
  initial_code <- if (name == "dlpfc") {
    source$initial[, scenario_index]
  } else {
    .colorectal_corruption(source, scenario_index)
  }
  sample_values <- if (name == "dlpfc") {
    source$samples
  } else {
    rep.int("CRC_Visium_HD", nrow(source$xy))
  }
  .validate_class_anchors(
    source$truth, initial_code, sample_values,
    context = paste0("Bundled ", if (name == "crc") "CRC" else name,
                     " benchmark scenario")
  )
  labels <- factor(
    source$levels[initial_code], levels = source$levels
  )
  truth <- factor(source$levels[source$truth], levels = source$levels)
  metadata <- c(
    source$metadata,
    as.list(scenario_design[1L, , drop = FALSE])
  )
  benchmark <- spatial_benchmark(
    xy = source$xy,
    labels = labels,
    truth = truth,
    samples = sample_values,
    boundary = source$boundary,
    regions = truth,
    sparse = source$sparse,
    name = as.character(scenario_design$scenario_id),
    metadata = metadata
  )
  benchmark$spot_id <- source$spot_id
  if (name == "dlpfc") {
    benchmark$subject <- source$subject
    benchmark$nearest_adjacent_layer <- source$nearest_adjacent_layer
  }
  benchmark$scenario <- scenario_design
  benchmark
}

.colorectal_corruption <- function(source, scenario_index) {
  condition <- source$design[scenario_index, , drop = FALSE]
  n <- nrow(source$xy)
  classes <- length(source$levels)
  count <- min(n - 1L, as.integer(round(condition$noise * n)))

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(1040000L + scenario_index)
  anchor_proximity <- as.numeric(!source$boundary)

  if (condition$mechanism == "random") {
    initial_rows <- sample.int(n, count)
    rows <- .select_class_preserving_corruption(
      c(initial_rows, setdiff(seq_len(n), initial_rows)), count, source$truth,
      boundary_proximity = anchor_proximity
    )
    replacement <- sample.int(classes - 1L, count, replace = TRUE)
    replacement <- replacement + (replacement >= source$truth[rows])
  } else if (condition$mechanism == "boundary") {
    candidates <- which(source$boundary)
    initial_rows <- if (length(candidates) >= count) {
      sample(candidates, count)
    } else {
      c(candidates, sample(setdiff(seq_len(n), candidates), count - length(candidates)))
    }
    rows <- .select_class_preserving_corruption(
      c(initial_rows, setdiff(seq_len(n), initial_rows)), count, source$truth,
      boundary_proximity = anchor_proximity
    )
    replacement <- source$adjacent[rows]
  } else {
    unit_xy <- sweep(source$xy, 2L, apply(source$xy, 2L, min), "-")
    unit_xy <- sweep(
      unit_xy, 2L, pmax(apply(unit_xy, 2L, max), .Machine$double.eps), "/"
    )
    if (condition$mechanism == "patch") {
      centers <- sample.int(n, 16L)
      distances <- vapply(
        centers, function(center) rowSums((unit_xy - unit_xy[center, ])^2),
        numeric(n)
      )
      owner <- max.col(-distances, ties.method = "first")
      distance <- distances[cbind(seq_len(n), owner)]
      priority <- order(distance)
      rows <- .select_class_preserving_corruption(
        priority, count, source$truth, boundary_proximity = anchor_proximity
      )
      replacement <- source$adjacent[centers][owner[rows]]
      same <- replacement == source$truth[rows]
      replacement[same] <- source$adjacent[rows[same]]
    } else {
      center <- sample.int(n, 1L)
      distance <- rowSums((unit_xy - unit_xy[center, ])^2)
      priority <- order(distance)
      rows <- .select_class_preserving_corruption(
        priority, count, source$truth, boundary_proximity = anchor_proximity
      )
      replacement <- rep.int(source$adjacent[center], count)
      same <- replacement == source$truth[rows]
      replacement[same] <- source$adjacent[rows[same]]
    }
  }

  initial <- source$truth
  initial[rows] <- replacement
  .validate_class_anchors(
    source$truth, initial, context = "CRC benchmark corruption"
  )
  initial
}

#' Human DLPFC Visium refinement benchmark
#'
#' A compact redistribution of spatial coordinates and manually annotated
#' cortical layers for 47,329 spots from 12 sections and three donors in the
#' Maynard et al. human DLPFC study. It also contains 45 frozen, anatomically
#' constrained label-corruption scenarios. Expression counts and tissue images
#' are deliberately excluded.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{xy}{A 47,329 by 2 coordinate matrix.}
#'   \item{truth}{Integer reference-layer codes.}
#'   \item{initial}{A 47,329 by 45 matrix of corrupted label codes.}
#'   \item{levels}{Seven cortical-layer names.}
#'   \item{samples}{Section identifiers.}
#'   \item{subject}{Donor identifiers.}
#'   \item{boundary}{Logical boundary indicator.}
#'   \item{sparse}{Logical sparse-layer indicator.}
#'   \item{design}{Scenario corruption design.}
#'   \item{metadata}{Source, license, and citation information.}
#' }
#' @source The `spatialLIBD` Bioconductor data package,
#'   \url{https://bioconductor.org/packages/spatialLIBD}.
#' @references
#' Maynard KR et al. (2021). Transcriptome-scale spatial gene expression in the
#' human dorsolateral prefrontal cortex. Nature Neuroscience 24, 425-436.
#'
#' Pardo B et al. (2022). spatialLIBD: an R/Bioconductor package to visualize
#' spatially-resolved transcriptomics data. BMC Genomics 23, 434.
"dlpfc_benchmark"

#' Human CRC Visium HD refinement benchmark
#'
#' A compact derivative of the 10x Genomics human colorectal cancer (CRC) Visium
#' HD dataset containing 194,541 annotated locations, two-dimensional coordinates,
#' 19 author-derived WSI region labels, and 60 deterministic corruption recipes.
#' Expression counts and tissue imagery are excluded. The derivative records all
#' transformations and is redistributed under CC BY 4.0.
#'
#' @format A list with coordinates, integer reference labels, adjacent-region
#'   labels, boundary and sparse indicators, corruption design, and attribution
#'   metadata.
#' @source 10x Genomics, Visium HD Spatial Gene Expression Library, Human
#'   Colorectal Cancer (FFPE), analyzed with Space Ranger 4.0.1,
#'   \url{https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc-v4}.
#' @references 10x Genomics (2025). Visium HD Spatial Gene Expression Library,
#'   Human Colorectal Cancer (FFPE). Published July 3, 2025. CC BY 4.0.
"colorectal_benchmark"
