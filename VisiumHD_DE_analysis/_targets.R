# Pipeline entry point for the targets workflow.
# Run `targets::tar_make()` to execute or resume the pipeline.
# See README.md for a full description of the analysis.
library(targets)
library(tarchetypes)

# ---------------------------------------------------------------------------
# Crew controller setup
# ---------------------------------------------------------------------------
# Heavy targets are dispatched to SLURM; lightweight bookkeeping runs locally.
# Controller names encode resource requests: slurm_{cpus}c{memory_gb}g.
# Assign a controller to a target via:
#   resources = tar_resources(crew = tar_resources_crew(controller = "slurm_1c40g"))
# ---------------------------------------------------------------------------
controller_local <- crew::crew_controller_local(
  name = "controller_local",
  workers = 3,
  retry_tasks = FALSE,
  seconds_idle = 10
)

create_slurm_controller <- function(
    cpus, memory_gb, name, workers = 10,
    seconds_idle = 10, seconds_timeout = 60,
    log_prefix = "crew_log",
    script_lines = "module load R/flexiblas/4.4.1 curl") {
  if (missing(name)) {
    name <- sprintf("slurm_%dc%dg", cpus, memory_gb)
  }
  crew.cluster::crew_controller_slurm(
    name = name,
    workers = workers,
    seconds_idle = seconds_idle,
    seconds_timeout = seconds_timeout,
    # retry_tasks = FALSE,
    options_cluster = crew.cluster::crew_options_slurm(
      script_lines = script_lines,
      memory_gigabytes_required = memory_gb,
      cpus_per_task = cpus,
      log_output = file.path("logs", sprintf("%s_%%A.txt", log_prefix)),
      log_error = file.path("logs",sprintf("%s_%%A.txt", log_prefix))
    )
  )
}

controllers <- do.call(crew::crew_controller_group,
  list(
    controller_local,
    create_slurm_controller(1, 20),
    create_slurm_controller(2, 40),
    create_slurm_controller(1, 40),
    create_slurm_controller(1, 80),
    create_slurm_controller(4, 80),
    create_slurm_controller(1, 120),
    create_slurm_controller(4, 120),
    create_slurm_controller(32, 500),
    create_slurm_controller(1, 500),
    create_slurm_controller(1, 1000)
  )
)

# ---------------------------------------------------------------------------
# Global pipeline options
# ---------------------------------------------------------------------------
# packages: loaded on every worker before executing a target.
# format = "qs": serialize results with the qs package (fast + compact).
# storage/retrieval = "worker": workers write/read their own results directly,
#   avoiding shipping large objects back through the main process.
# ---------------------------------------------------------------------------
tar_option_set(
  packages = c(
    "tibble", "SpatialExperiment", "RColorBrewer", "gridExtra",
    "grDevices", "SingleCellExperiment", "ggplot2", "tidyverse",
    "scater", "scran", "ggspavis", "limma", "edgeR", "cowplot",
    "patchwork", "tweeDEseqCountData"
  ),
  format = "qs",
  # memory = "transient", garbage_collection = TRUE,
  storage = "worker", retrieval = "worker",
  error = "stop",
  controller = controllers,
  resources = tar_resources(
    crew = tar_resources_crew(controller = "controller_local"),
  )
)
# Source all helper functions in R/ into the pipeline environment.
tar_source()

# Auto-discover and evaluate every *.R file in targets/.
# Each file must return a list() of tar_target() calls; those lists are
# concatenated here to form the complete set of pipeline targets.
unlist(
  lapply(
    X = list.files(
      path = file.path("targets"),
      pattern = "\\.R$",
      full.names = TRUE,
      all.files = TRUE,
      recursive = TRUE
    ),
    FUN = function(path) {
      eval(
        expr = parse(text = readLines(con = path, warn = FALSE)),
        envir = targets::tar_option_get(name = "envir")
      )
    }
  )
)
