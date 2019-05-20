# =============================================================================
# 
# Dynamic wildlife occupancy models using automated acoustic monitoring data
# C. Balantic & T. Donovan
# Submitted to Ecological Applications December 2018
# Script for running simulation experiment and bias calculations
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
#       SUPPLEMENTAL DATA, AND WILL LIKELY TAKE A FEW DAYS TO RUN, DEPENDING 
#       ON YOUR MACHINE. Navigate to SECTION 3 to reproduce visualizations using 
#       the supplemental data provided. 

exp.folder <- 'Data/simulation-results/' # experiment folder
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
f.pct.confirms <- c(0.025, 0.05) # manual confirmation effort levels
f.pr.TP <- c('high_low', 'med_med') # classifier types: good or bad

# flag for whether to eliminate any detections below 0.5 probability
# ==> this is the equivalent of applying the classification principle in chapter 1
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
  rep.scenario = 500                # n replicates for each scenario
)


# CALCULATE BIAS USING RPRESENCE USING THE MILLER 2013 MODEL:
amdata <- readRDS(paste0(exp.folder, 'amdata_object.RDS'))
bias.return <- simCalcBias(amdata = amdata,
                           psi.formula = 'psi ~ 1', 
                           gamma.formula = 'gamma ~ 1',
                           epsilon.formula = 'epsilon ~ 1',
                           p11.formula = 'p11 ~ 1',
                           p10.formula = 'p10 ~ 1',
                           b.formula = 'b ~ 1',
                           model.name = 'intercept',
                           experiment.folder = exp.folder,
                           save.detections = TRUE,
                           use.informed.inits = TRUE) 




# SECTION 3: VISUALIZING RESULTS ==============================================

exp.folder <- 'Data/simulation-results/' # experiment folder

# LOOK AT STATE PARAMETERS (a la Fig. 6 & Appendix S2 Fig. S1)
simPlotBias(
  experiment.folder = exp.folder,
  false.alarm.rate.paper = 0.8,
  pct.confirms.paper = c(0.025, 0.05),
  n.days.paper = c(1, 3),
  any.true.thresholds.paper = c(0.8, 0.95),
  plot.state.parameters = TRUE,
  remove.warnings = 'convergence',
  zoom = FALSE,
  scale.y.continuous = scale_y_continuous(
    breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
    limits = c(-0.5, 0.5))
)


# STATE PARAMETERS ZOOMED IN TO MOST CONSERVATIVE CONDITIONS (a la Fig. 7)
simPlotBias(
  experiment.folder = exp.folder,
  false.alarm.rate.paper = 0.8,
  pct.confirms.paper = c(0.025, 0.05),
  n.days.paper = 1, # only look at most conservative condition
  any.true.thresholds.paper = 0.95, # only look at most conservative condition
  plot.state.parameters = TRUE,
  remove.warnings = 'convergence',
  zoom = TRUE,
  scale.y.continuous = scale_y_continuous(
    breaks = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15),
    limits = c(-0.15, 0.15))
)


# LOOK AT DETECTION PARAMETERS (a la Fig. 8 & Appendix S2 Fig. S2)
simPlotBias(
  experiment.folder = exp.folder,
  false.alarm.rate.paper = 0.8,
  pct.confirms.paper = c(0.025, 0.05),
  n.days.paper = c(1, 3),
  any.true.thresholds.paper = c(0.8, 0.95),
  plot.state.parameters = FALSE,
  remove.warnings = 'convergence',
  zoom = FALSE,
  scale.y.continuous = scale_y_continuous(
    breaks = c(-0.05, -0.02, 0, 0.02, 0.05),
    limits = c(-0.05, 0.05))
)



