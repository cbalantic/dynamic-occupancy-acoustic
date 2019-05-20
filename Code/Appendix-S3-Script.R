# =============================================================================
# 
# Dynamic wildlife occupancy models using automated acoustic monitoring data
# C. Balantic & T. Donovan
# Submitted to Ecological Applications December 2018
# APPENDIX S3: Classic Mackenzie 2003 Dynamic Model Results
# Script for running appendix experiment and bias calculations
# for the Mackenzie 2003 dynamic occupancy model
#
# =============================================================================

# SECTION 1: SET UP ===========================================================
library(data.table)
library(AMModels)
library(RPresence)
library(ggplot2)
library(gridExtra)
library(grid)

# Download the "Code" and "Data" folders, then set your working directory appropriately.

# Source in all functions from the simulation-functions.R file
source('Code/Simulation-Functions.R')

# Install PRESENCE software and R presence here: 
# https://www.mbr-pwrc.usgs.gov/software/presence.html
# Simulations were run with RPresence version 12.10 and R version 3.4.4

# SECTION 2: WRAPPER FUNCTION FOR SIMULATION + BIAS CALCULATIONS ==============
# NOTE: RUNNING CODE IN THIS SECTION WILL OVERWRITE EXISTING PROVIDED 
#       SUPPLEMENTAL DATA, AND WILL LIKELY TAKE SEVERAL HOURS TO RUN, DEPENDING 
#       ON YOUR MACHINE. Navigate to SECTION 3 to reproduce visualizations using 
#       the supplemental data provided.  

# Set up the simulation parameters:
exp.folder <- 'Data/appendix-results/' # experiment folder
n.seasons <- 2 # simulate for two seasons
n.days.in.season <- 30 # 30 days in each season
occ.probs <- c(0.2, 0.6) # probabilities of occupancy
col.probs <- c(0.05, 0.25) # probabilities of colonization
ext.probs <- c(0.05, 0.25) # probabilities of extinction
voc.rates <- c(20, 100) # hourly rates, songs per hour.
f.n.sites <- 100 # number of sites
f.daily.n.samples <- 5 # number of recordings to take each day
f.false.alarm.rates <- 0.8 # underlying false alarm rate in the soundscape
f.any.true.thresholds <- c(0.8, 0.95) # survey-level detection thresholds
f.n.days.in.survey <- c(1, 3) # survey aggregation lengths (1 or 3 days)
f.pct.confirms <- 0 # 0% means we use the Mackenzie 2003 model
f.pr.TP <- c('high_low', 'med_med') # classifier types: good or bad

# flag for whether to eliminate any detections below 0.5 probability
f.cull.FP <- FALSE

# CREATE SCENARIO GRID:
s.grid <- simScenarioGrid(
  
  # Models:
  occupancy.probs = occ.probs,
  colonization.probs = col.probs,
  extinction.probs = ext.probs,
  vocalization.rates = voc.rates,
  
  # Season parameters:
  n.seasons = n.seasons,
  n.days.in.season = n.days.in.season,
  
  # ==> Dynamics Dials:
  n.sites = f.n.sites,
  
  # ==> Miller 2013 dials
  pct.confirm = f.pct.confirms,
  
  # ==> Aggregation Dials:
  any.true.thresholds = f.any.true.thresholds,
  n.days.in.survey = f.n.days.in.survey,
  daily.n.samples = f.daily.n.samples,
  false.alarm.rates = f.false.alarm.rates,
  pr.TP.shapes = f.pr.TP,
  cull.FP = f.cull.FP
)


# RUN SIMULATIONS FOR ALL EXPERIMENTS
set.seed(3)
simOccupancyExperiment(
  scenario.grid = s.grid,           # input the grid of scenarios to be run
  experiment.folder = exp.folder,   # folder for results storage
  voc.time.units = 'hourly',        # model lambda rates on an hourly basis
  rep.scenario = 100                # n replicates for each scenario: only 100 reps for appendix expt
)


# CALCULATE BIAS USING RPRESENCE USING THE MILLER 2013 MODEL:
amdata <- readRDS(paste0(exp.folder, 'amdata_object.RDS'))
bias.return <- simCalcBias(amdata = amdata,
                           psi.formula = 'psi ~ 1', 
                           gamma.formula = 'gamma ~ 1',
                           epsilon.formula = 'epsilon ~ 1',
                           p.formula = 'p ~ 1',
                           model.name = 'intercept',
                           experiment.folder = exp.folder,
                           save.detections = TRUE,
                           use.informed.inits = TRUE) 

# SECTION 3: VISUALIZING RESULTS ==============================================

exp.folder <- 'Data/appendix-results/' # experiment folder

# LOOK AT STATE PARAMETERS
simPlotBias(
  experiment.folder = exp.folder,
  false.alarm.rate.paper = 0.8,
  pct.confirms.paper = 0,                   # Look only at 0% confirmation (Mackenzie 2003 model)
  n.days.paper = c(1,3),                    # Look at both aggregations
  any.true.thresholds.paper = c(0.8, 0.95), # Look at both thresholds
  plot.state.parameters = TRUE,             # Plot state parameters
  remove.warnings = 'convergence',
  zoom = FALSE,
  scale.y.continuous = scale_y_continuous(
    breaks = seq(from = -1, to = 1, by = 0.25),
    limits = c(-1, 1))
)

# LOOK AT DETECTION PARAMETERS
simPlotBias(
  experiment.folder = exp.folder,
  false.alarm.rate.paper = 0.8,
  pct.confirms.paper = 0,                   # Look only at 0% confirmation (Mackenzie 2003 model)
  n.days.paper = c(1,3),                    # Look at both aggregations
  any.true.thresholds.paper = c(0.8, 0.95), # Look at both thresholds
  plot.state.parameters = FALSE,            # Plot detection parameters
  remove.warnings = 'convergence',
  zoom = FALSE,
  scale.y.continuous = scale_y_continuous(
    breaks = seq(from = -1, to = 1, by = 0.25),
    limits = c(-1, 1))
)

