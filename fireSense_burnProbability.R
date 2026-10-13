# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_burnProbability",
  description = paste(
    "Estimates each pixel's probability of burning in a given year by Monte Carlo:",
    "repeatedly re-draws ignitions, escapes, and fire spread (via",
    "fireSense_ignitionPredict and fireSense_burn) without touching the",
    "simulation's realized fire outputs, and tallies how often each pixel burns."),
  keywords = c("fire", "probability", "monte carlo", "fireSense"),
  authors = c(
    person("Dominique", "Caron", email = "", role = c("aut", "cre"))
  ),
  childModules = character(),
  version = list(fireSense_burnProbability = "0.1.0", SpaDES.core = "0.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "fireSense_burnProbability.Rmd"),
  reqdPkgs = list("data.table", "terra", "SpaDES.core (>=3.0.4)"),
  loadOrder = list(after = c("fireSense_ignitionPredict", "fireSense_spreadPredict", "fireSense_burn")),
  parameters = bindrows(
    defineParameter("nIterations", "integer", 10000L, 1L, NA,
                    desc = paste("Number of Monte Carlo iterations (independent ignition/escape/spread",
                                 "draws) used to estimate burn probability at each `.runInterval`.")),
    defineParameter(".runInitialTime", "numeric", default = start(sim), NA, NA,
                    desc = "Time of the first burn-probability estimate."),
    defineParameter(".runInterval", "numeric", default = 1, NA, NA,
                    desc = paste("Years between burn-probability estimates. `NA` estimates once, at",
                                 "`.runInitialTime`. Each estimate uses that year's ignition/spread",
                                 "covariates and model predictions, held fixed across all `nIterations` draws.")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used.")
  ),
  inputObjects = bindrows(
    expectsInput("flammableRTM", "SpatRaster", sourceURL = NA,
                 "Binary raster, 1 where the pixel is flammable. Template for `burnProbability`."),
    expectsInput("rasterToMatch", "SpatRaster", sourceURL = NA,
                 "Template raster for the study area (as used by `fireSense_burn`)."),
    expectsInput("fireSense_SpreadPredicted", "SpatRaster",
                 "This year's spread probability surface, from `fireSense_spreadPredict`."),
    expectsInput("fireSense_SpreadSD", "SpatRaster|numeric",
                 "This year's fitted year-random-effect sd on logit spread probability, from `fireSense_spreadPredict`."),
    expectsInput("ignitionsAndEscapes", "data.table",
                 paste("This year's *realized* ignitions/escapes, from `fireSense_ignitionPredict`.",
                       "Snapshotted and restored; never left modified by the Monte Carlo loop.")),
    expectsInput("nonEscapedFireSizesHa", "numeric",
                 "Sizes (ha) of observed small fires, from `fireSense_dataPrepFit`; as used by `fireSense_burn`."),
    expectsInput("burnMap", "SpatRaster",
                 paste("This year's *realized* cumulative burn map, from `fireSense_burn`.",
                       "Snapshotted and restored, never left modified by the Monte Carlo loop.")),
    expectsInput("rstCurrentBurn", "SpatRaster",
                 "This year's *realized* current-year burn indicator, from `fireSense_burn`. Snapshotted and restored."),
    expectsInput("rstAnnualBurnID", "SpatRaster",
                 "This year's *realized* annual fire IDs, from `fireSense_burn`. Snapshotted and restored."),
    expectsInput("burnDT", "data.table",
                 "This year's *realized* burned-pixel table, from `fireSense_burn`. Snapshotted and restored."),
    expectsInput("burnSummary", "data.table",
                 paste("Cumulative per-fire summary, from `fireSense_burn`. Snapshotted and restored",
                       "(each Monte Carlo draw would otherwise append a spurious row)."))
  ),
  outputObjects = bindrows(
    createsOutput("burnProbability", "SpatRaster",
                  desc = paste("Estimated probability each flammable pixel burns this interval:",
                               "fraction of `nIterations` Monte Carlo draws in which the pixel burned.",
                               "`NA` where not flammable."))
  )
))

#' Event dispatcher
#'
#' Events: `init` (schedules the first `estimate`), `estimate` (runs the Monte
#' Carlo loop and reschedules itself every `.runInterval`).
#'
#' @param sim A `simList`.
#' @param eventTime Time of the event.
#' @param eventType Name of the event.
#' @param debug Not used.
#'
#' @return The `simList`, invisibly.
doEvent.fireSense_burnProbability <- function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {
      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "estimate",
                           eventPriority = 5.14) ## after fireSense_burn (5.13)
    },
    estimate = {
      sim <- estimateBurnProbability(sim)

      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "estimate",
                             eventPriority = 5.14)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'",
                  sep = ""
    ))
  )

  invisible(sim)
}

#' Estimate this interval's burn probability by Monte Carlo
#'
#' Snapshots every `sim$` object that `fireSense_ignitionPredict::IgnitionPredictRun()`
#' and `fireSense_burn::burn()` write to, then runs `P(sim)$nIterations` independent
#' draws of ignitions/escapes and fire spread (re-calling those modules' own functions,
#' via `sim$.mods`, so the Monte Carlo draws use exactly the same stochastic logic as a
#' real simulation year), tallying how many iterations burn each pixel. The snapshotted
#' objects are restored before returning, so this module never alters the "realized"
#' simulation outputs other modules (e.g. carbon accounting) depend on -- it only adds
#' `burnProbability`.
#'
#' Deterministic spread probability (`fireSense_SpreadPredicted`) is *not* recomputed
#' each iteration: it does not change between draws within one interval, only the
#' ignition/escape draws and the year random effect / spread realization do.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly, with `burnProbability`.
estimateBurnProbability <- function(sim) {
  nIter <- P(sim)$nIterations

  snapshotNames <- c("ignitionsAndEscapes", "burnMap", "rstCurrentBurn", "rstAnnualBurnID",
                     "burnDT", "burnSummary")
  snapshot <- snapshotSimObjects(sim, snapshotNames)

  burnCount <- rast(sim$flammableRTM)
  burnCount[] <- 0
  flamVals <- as.vector(sim$flammableRTM)
  burnCount[flamVals == 0] <- NA

  for (i in seq_len(nIter)) {
    sim <- runAsModule(sim, "fireSense_ignitionPredict", sim$.mods$fireSense_ignitionPredict$IgnitionPredictRun)
    sim <- runAsModule(sim, "fireSense_burn", sim$.mods$fireSense_burn$burn)

    burned <- !is.na(terra::values(sim$rstCurrentBurn, mat = FALSE)) &
      terra::values(sim$rstCurrentBurn, mat = FALSE) == 1
    burnCount[which(burned)] <- burnCount[which(burned)] + 1
  }

  sim <- restoreSimObjects(sim, snapshot)

  sim$burnProbability <- burnCount / nIter

  invisible(sim)
}

#' Run a delegate module's function with `current(sim)` set to that module
#'
#' `P(sim)` (no explicit module name) resolves against `sim@current$moduleName`,
#' i.e. whichever module the scheduler currently has active -- not whichever
#' module's function is lexically executing. Calling e.g.
#' `sim$.mods$fireSense_burn$burn(sim)` directly from inside
#' `fireSense_burnProbability`'s own event leaves `sim@current$moduleName` set
#' to `"fireSense_burnProbability"`, so any `P(sim)$<param>` lookups inside
#' `burn()` would silently return `NULL` instead of `fireSense_burn`'s
#' parameters. This wraps such a call so `sim@current$moduleName` is
#' temporarily set to `moduleName`, restoring the prior value on the
#' `simList` that `fn` returns.
#'
#' @param sim A `simList`.
#' @param moduleName character; the name of the module whose function is being run.
#' @param fn The module's event function (e.g. `sim$.mods$fireSense_burn$burn`).
#'
#' @return The `simList`, as returned by `fn`.
runAsModule <- function(sim, moduleName, fn) {
  prevCurrent <- sim@current
  sim@current <- modifyList(prevCurrent, list(moduleName = moduleName))
  sim <- fn(sim)
  sim@current <- prevCurrent
  sim
}

#' Snapshot a set of `sim$` objects before they are overwritten
#'
#' `SpatRaster`s are deep-copied (`terra::deepcopy()`): otherwise they alias the same
#' underlying pointer, and `restoreSimObjects()` would "restore" a reference to data
#' the Monte Carlo loop has already overwritten. `data.table`s are copied with
#' `data.table::copy()` for the same reason (they are modified by reference).
#'
#' @param sim A `simList`.
#' @param names character; the `sim$` object names to snapshot.
#'
#' @return A named list holding each object's copy (or `NULL` if absent).
snapshotSimObjects <- function(sim, names) {
  setNames(lapply(names, function(n) {
    obj <- sim[[n]]
    if (is.null(obj)) return(NULL)
    if (inherits(obj, "SpatRaster")) return(terra::deepcopy(obj))
    if (is.data.table(obj)) return(data.table::copy(obj))
    obj
  }), names)
}

#' Restore `sim$` objects from a snapshot
#'
#' @param sim A `simList`.
#' @param snapshot A named list, as returned by `snapshotSimObjects()`.
#'
#' @return The `simList`, invisibly, with each snapshotted object restored.
restoreSimObjects <- function(sim, snapshot) {
  for (n in names(snapshot)) sim[[n]] <- snapshot[[n]]
  invisible(sim)
}

#' Supply default inputs
#'
#' There are none.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly.
.inputObjects <- function(sim) {
  return(invisible(sim))
}
