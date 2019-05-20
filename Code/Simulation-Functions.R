# =============================================================================
# 
# Dynamic wildlife occupancy models using automated acoustic monitoring data
# C. Balantic & T. Donovan
# Submitted to Ecological Applications December 2018
# Support functions for running simulation and bias calculations
#
# =============================================================================

# simAggregate ================================================================
#' Aggregate multiple recording sessions into surveys to produce encounter histories for use in RPresence.
#' @name simAggregate
#' @title Aggregate multiple recording sessions into surveys to produce encounter histories for use in RPresence
#' @param dynamics 'dynamics' object output from simDynamics()
#' @param vocalizations A list of vocalizations produced by simVocals()
#' @param detections A list of detections produced by simDetections()
#' @param pct.confirm Percentage of total surveys to confirm. For now, confirmation surveys are randomly allocated across sites and seasons.
#' @param n.days.in.survey Number of days that should constitute one survey
#' @param any.true.threshold Threshold of the probability of any detection being true above which to designate a 1 (detection) for the survey period
#' @param cull.FP Logical flag for whether to cull from the dataset any detections with target signal probability (prTP) below 0.5.
#' @return A list containing two data.tables, one of encounter histories, and one of aggregated probabilities.
#' @details TBD.
#' @seealso Coming soon... link to simDynamics, simVocalizations, simDetections
#' @export
#' @examples
#'
#' \dontrun{
#'
#' }
#'

simAggregate <- function(dynamics,
                         vocalizations,
                         detections,
                         pct.confirm,
                         n.days.in.survey = 1,
                         any.true.threshold = 0.99,
                         cull.FP = FALSE
)
{
  
  # Note: currently, function does not do anything to reconcile uneven numbers of days into survey groups.
  
  # Set up parameters
  n.days.in.survey <- as.numeric(n.days.in.survey)
  any.true.threshold <- as.numeric(any.true.threshold)
  n.seasons <- length(grep(x = names(dynamics), pattern = 'season'))
  n.sites <- length(unique(dynamics$site))
  n.days.in.season <- length(unique(vocalizations$season.day))
  daily.n.samples <- length(unique(vocalizations$recording))
  surveys.per.season <- floor(n.days.in.season/n.days.in.survey)
  
  # Initialize matrix to hold encounter histories
  unconfirmed.EH <- matrix(data = 0, nrow = n.sites, ncol = surveys.per.season*n.seasons)
  colnames(unconfirmed.EH) <- unlist(lapply(1:n.seasons, function(x)
    paste0(x, '-', 1:surveys.per.season)))
  
  # Assign indices of which surveys to confirm:
  confirm.inds <- sort(sample(x = length(unconfirmed.EH),
                              size = length(unconfirmed.EH)*pct.confirm, replace = FALSE))
  to.confirm.EH <- confirmed.EH <- unconfirmed.EH
  to.confirm.EH[confirm.inds] <- 1
  
  # Anything unconfirmed will be -1. 1 will be a TP, 0 will be an FP
  confirmed.EH[confirmed.EH == 0] <- -1
  
  # For each season
  seasons.ag <- list()
  for (seas in 1:n.seasons) {
    
    # Find columns that indicate seasons
    season.cols <- grep(x = colnames(unconfirmed.EH), pattern = paste0(seas,'-'),
                        fixed = TRUE)
    
    sites.ag <- list()
    
    # for each location
    for (L in 1:n.sites) {
      
      # Read in detections for this site and season
      dets.site.season <- detections[site == L & season == seas]
      
      # Cull events with target signal probability below 0.5, if desired
      if (cull.FP == TRUE) {
        dets.site.season <- dets.site.season[pr.TP >= 0.5]
      }
      
      # Read in vocalization data for this site and season:
      vocs.site.season <- vocalizations[site == L & season == seas]
      
      # Create survey groups for each sampled day:
      survey.groups <- ceiling(vocs.site.season$season.day/n.days.in.survey)
      
      # If uneven n.days.in.survey, make corrections and tack extra surveys into the last group
      survey.groups[survey.groups > surveys.per.season] <-
        surveys.per.season
      vocs.site.season[, survey := survey.groups]
      
      # If any detections occurred at this site for this season:
      if (nrow(dets.site.season) > 0) {
        
        # Calculate the probability that a detection is a false alarm
        dets.site.season[, pr.FP := 1 - pr.TP]
        
        # Merge detections and vocalizations
        dets.site.season <- merge(x = dets.site.season,
                                  y = vocs.site.season[,c('season.day', 'recording', 'survey')],
                                  by = c('season.day', 'recording'),
                                  all.x = TRUE)
        
        # Aggregate false alarms and target signals, calculate the 
        # probably that there is ANY target signal within the survey period
        agFP <- dets.site.season[ , prod(pr.FP, na.rm = TRUE), by = survey]
        agTP <- dets.site.season[ , prod(pr.TP, na.rm = TRUE), by = survey]
        ag <- agFP[agTP, on = 'survey'] 
        colnames(ag) <- c('survey', 'p.all.false', 'p.all.true')
        ag[ , p.any.true := 1 - p.all.false]
        ag[ , p.any.false := 1 - p.all.true]
        ag[ , unconfirmed.EH := as.integer(p.any.true > any.true.threshold)]
        # p.all.false ==> product of prFP across detected events for each survey
        # p.all.true ==> product of prTP across detected events for each survey
        # p.any.true ==> 1 - p.all.false
        # p.any.false ==> 1 - p.all.true
        
        # Figure out if any surveys were missed AT THE RECORDING LEVEL
        # (not at the detection level); assign NA if so.
        missing.surveys <- which(!(1:surveys.per.season %in%
                                     unique(vocs.site.season$survey)))
        if (length(missing.surveys) > 0) {
          missing.surveys.dt <- data.table(survey = missing.surveys,
                                           p.all.false = as.numeric(NA),
                                           p.all.true = as.numeric(NA),
                                           p.any.true = as.numeric(NA),
                                           p.any.false = as.numeric(NA),
                                           unconfirmed.EH = as.integer(NA)) # this is NA because the survey was missed completely
          ag <- rbind(ag, missing.surveys.dt)
          setkey(ag, survey)
        }
        
        # If surveys occurred but resulted in no detections, add these to enc history:
        survs.w.no.detection <- which(!(1:surveys.per.season %in% ag$survey))
        if (length(survs.w.no.detection) > 0) {
          no.dets.dt <- data.table(survey = survs.w.no.detection,
                                   p.all.false = as.numeric(NA),
                                   p.all.true = as.numeric(NA),
                                   p.any.true = as.numeric(NA),
                                   p.any.false = as.numeric(NA),
                                   unconfirmed.EH = as.integer(0)) # 0 because nothing detected
          ag <- rbind(ag, no.dets.dt)
          setkey(ag, survey)
        }
      } # end if detections occurred
      
      # If no detections occurred at this site for this season,
      #  record a bunch of 0s!:
      if (nrow(dets.site.season) == 0) {
        ag <- data.table(survey = 1:surveys.per.season,
                         p.all.false = as.numeric(NA),
                         p.all.true = as.numeric(NA),
                         p.any.true = as.numeric(NA),
                         p.any.false = as.numeric(NA),
                         unconfirmed.EH = as.integer(0)) # 0 because nothing detected
      }
      
      # Store the results
      sites.ag[[L]] <- ag
      unconfirmed.EH[L, season.cols] <- ag$unconfirmed.EH
      
      
      # Gather detection histories for multiple detection state specification
      confirm.survs <- colnames(to.confirm.EH[L, season.cols, drop = FALSE])[
        which(to.confirm.EH[L, season.cols] == 1)]
      
      if (length(confirm.survs) > 0) {
        cs <- as.numeric(
          unlist(
            lapply(
              strsplit(x = confirm.survs, split = '-', fixed = TRUE),
              '[[',
              2)))
        
        # VERIFY the surveys that should be confirmed (if there are any detections)
        # if no detections, site-survey will remain as 0 and count as checked (there wont be any FP or TP to confirm)
        if (nrow(dets.site.season) > 0) {
          dv <- merge(x = vocs.site.season[,c('season.day', 'recording', 'survey')],
                      y = dets.site.season,
                      by = c('season.day', 'recording', 'survey'),
                      all.x = TRUE)
          d <- dv[survey %in% cs]
          
          
          # there may be NAs -- that means we were recording but there was no detection
          # if we have NAs, replace them with false.alarm
          # to show that we failed to detect anything for this survey
          # (this way no "true" gets added for the encounter history)
          d[is.na(vocal.type), vocal.type := 'false.alarm']
          assign2 <- d[,sum(vocal.type == 'true.vocal'),
                       by = survey][,V1 > 0]
          
          names(assign2) <- confirm.survs
          confirmed.EH[L, which(colnames(confirmed.EH) %in% confirm.survs)] <- assign2
        }
        
      }# end if confirm.survs > 0
    } # end for location
    
    # seasons.ag[[seas]] <- sites.ag # dont' need to save this
    
  } # end for season
  
  # Reconcile the uncertain detections (unconfirmed.EH) with the certain detections (confirmed.EH)
  EH <- unconfirmed.EH
  EH[confirmed.EH == 1] <- 2
  EH[confirmed.EH == 0] <- 0
  EH
  
}

# simCalcBias =================================================================
#' Using RPresence::occMod(), fit 'do.fp' models from simOccupancyExperiment() outputs, and compare estimates to simulated truth.
#' @name simCalcBias
#' @title Using RPresence::occMod(), fit 'do.fp' models from simOccupancyExperiment() outputs, and compare estimates to simulated truth.
#' @param amdata Input amdata object produced by simOccupancyExperiment
#' @param psi.formula String model formula for psi (occupancy) estimate for occMod(type = 'do.fp') model fit. String must begin with 'psi ~'.
#' @param gamma.formula String model formula for gamma (colonization) estimate for occMod(type = 'do.fp') model fit. String must begin with 'gamma ~'.
#' @param epsilon.formula String model formula for epsilon (extinction) estimate for occMod(type = 'do.fp') model fit. String must begin with 'epsilon ~'.
#' @param p11.formula String model formula for p11 estimate for occMod(type = 'do.fp') model fit. String must begin with 'p11 ~'.
#' @param p10.formula String model formula for p10 estimate for occMod(type = 'do.fp') model fit. String must begin with 'p10 ~'.
#' @param b.formula String model formula for b estimate for occMod(type = 'do.fp') model fit. String must begin with 'b ~'.
#' @param p.formula String model formula for p estimate for occMod(type = 'do.1') model fit. String must begin with 'p ~'.
#' @param experiment.folder Path to store RDS results of each experiment (sends each exp to a folder separately for each scenario during troubleshooting stage... to at least store everything in case things fail. Temporary convenience for Cathleen in dev phase.)
#' @param model.name Character string of how to name the temporary RDS file (ie 'best', 'intercept')..temporary convenience for Cathleen.
#' @param save.detections Logical flag for whether to save AVERAGED bias of detection parameter estimates into a single table.
#' @param use.informed.inits Logical flag for whether to use informed inits. If yes, inputs inits into occMod() that match the simulated values in order to help convergence to a global minimum for the likelihood.
#' @param ... Additional arguments to pass to RPresence::occMod()
#' @return Data.table of bias estimates
#' @details tbd
#' @seealso Coming soon
#' @export
#' @examples
#'
#' #'\dontrun{
#'x <- NA
#'}
#

simCalcBias <- function(amdata,
                        psi.formula,
                        gamma.formula,
                        epsilon.formula,
                        p11.formula,
                        p10.formula,
                        b.formula,
                        p.formula,
                        experiment.folder,
                        model.name,
                        save.detections = FALSE,
                        use.informed.inits = FALSE,
                        # control = list(maxit = 2000), <=== e.g. how do i specify something like this for presence?
                        # method = 'BFGS'
                        ...)
{
  
  # Add a forward slash to experiment folder path if missing
  if (grepl("\\/$", experiment.folder) == FALSE) {
    experiment.folder <- paste0(experiment.folder, '/')
  }
  
  # Read in experiment data:
  scenario.grid <- amdata@data$scenario.grid
  survey.det <- amdata@data$detection.summary
  dyn <- amdata@data$dynamics
  enc.hist <- amdata@data$encounter.histories
  types <- c('occ', 'col', 'ext') # params being estimated
  
  # Loop through each scenario to calculate bias:
  for (ex in 1:length(dyn)) {
    
    cat('Scenario', ex, '\n')
    
    one.scen <- scenario.grid[ex]
    
    # Collect and save any convergence warnings
    # Establish a new warnings table if ex == 1
    if (ex == 1) {
      warns <- data.table(scenario = numeric(),
                          rep = numeric(),
                          conv.warning = character(),
                          VC.warning = logical())
    } else {
      # Otherwise, read in old warnings table and add to that one
      warns <- readRDS(file = paste0(experiment.folder,'presence_warnings.RDS'))
    }
    
    # If confirmation percentage is 0, we fit the regular dynamic occupancy model,
    #     which is FP-ignorant
    #    ('do.1' in PRESENCE lingo)
    # Otherwise, we fit the dynamic Miller model (with FPs)
    #    ('do.fp' in PRESENCE lingo)
    ifelse(one.scen$pct.confirm == 0, model.type <- 'do.1', model.type <- 'do.fp')
    
    # Gather helpful housekeeping info:
    scen <- one.scen$scenario       # current scenario
    one.dyn <- dyn[[ex]]            # these dynamics
    one.enc.hist <- enc.hist[[ex]]  # this encounter history
    n.seasons <- one.scen$n.seasons
    n.surveys.season <- ncol(one.enc.hist[[1]])/2
    
    # Calculate bias for each trial of this experiment epoch
    rep.scenario <- length(one.enc.hist)
    nrows <- rep.scenario*length(types)
    n.sites <- nrow(one.enc.hist[[1]])
    
    occ.est <- col.est <- ext.est <-
      matrix(0, nrow = n.sites,
             ncol = rep.scenario,
             dimnames = list(paste0('site.', 1:n.sites),
                             paste0('rep.', 1:rep.scenario)))
    
    # Gather the true simulated values --
    #  these can be pulled out of the i rep.scenario loop for speed since
    #  they'll be the same each time (because we are using an intercept model in the paper)
    #  (this will not work if a reader/reviewer tries to add any covariates to models)
    occ.actuals <- one.dyn[[1]]$pr.occ
    col.actuals <- one.dyn[[1]]$pr.col
    ext.actuals <- one.dyn[[1]]$pr.ext
    
    # Set up to collect detection parameter estimates
    #  p11, p10, and b only pertain to the dynamic miller model (do.fp)
    #  meanwhile, the regular dynamic model (do.1) only estimates 'p'
    p11.est <- p10.est <- b.est <- p.est <-
      matrix(0,
             nrow = n.sites*n.seasons*n.surveys.season,
             ncol = rep.scenario,
             dimnames = list(paste0('unit.', 1:(n.sites*n.seasons*n.surveys.season)),
                             paste0('rep.', 1:rep.scenario)))
    
    # These are different each rep.scenario, but we have already summarized them
    # in detection.summary, so we can pull them out right here as an average
    if (model.type == 'do.1') {
      p11.actuals <- survey.det[scenario == ex, p11] # these are also the p.actuals if using regular (non miller) model
    }
    if (model.type == 'do.fp') {
      p11.actuals <- survey.det[scenario == ex, p11] + survey.det[scenario == ex, b] # p11 is p + b
    }
    p10.actuals <- survey.det[scenario == ex, p10]
    b.actuals <- survey.det[scenario == ex, b]
    
    # For each replicate of this scenario, calculate bias
    for (i in 1:rep.scenario) {
      
      cat(i, '\n')
      
      single.rep.enc <- one.enc.hist[[i]]
      single.rep.dyn <- one.dyn[[i]]
      one.pao <- createPao(data = single.rep.enc,
                           nsurveyseason = rep(n.surveys.season, n.seasons))
      
      ifelse(model.type == 'do.fp',
             form.list <- list(psi.formula, gamma.formula, epsilon.formula,
                               p11.formula, p10.formula, b.formula),
             form.list <- list(psi.formula, gamma.formula, epsilon.formula,
                               p.formula))
      formulae <- lapply(form.list, as.formula)
      
      if (use.informed.inits) {
        # assumes no heterogeneity in simulation of params across sites
        #   e.g., same psi value was used to simulate occupancy of each site
        psi.init <- unique(occ.actuals)
        gamma.init <- unique(col.actuals)
        epsilon.init <- unique(ext.actuals)
        
        # If using the Mackenzie et al 2003 dynamic occupancy model:
        if (model.type == 'do.1') {
          # only need p11, aka p, detection param for do.1
          p.init <- mean(p11.actuals)
          
          ifelse(p.init == 1,
                 p.init <- p.init - 0.1, # bring within boundaries
                 p.init <- p.init)
          
          init.vals <- qlogis(c(psi.init, gamma.init, epsilon.init,
                                p.init))
        }
        
        # If using the Miller et al 2013 false positives dynamic occ model:
        if (model.type == 'do.fp') {
          # Detection parameters  need to be averaged instead
          #   bc there are difference across sites due to the nature of our sim
          p11.init <- mean(p11.actuals)
          p10.init <- mean(p10.actuals)
          b.init <- mean(b.actuals)
          
          # Check that detection param inits are not at a boundary condition (0|1)
          #  which will only happen in our simulation if using 100% confirmation
          #  if so, need to make them within boundaries
          ifelse(p11.init == 1,
                 p11.init <- p11.init - 0.1, # bring within boundaries
                 p11.init <- p11.init)
          
          ifelse(b.init == 1,
                 b.init <- b.init - 0.1, # bring within boundaries
                 b.init <- b.init)
          
          ifelse(p10.init == 0,
                 p10.init <- p10.init + 0.1, # bring within boundaries
                 p10.init <- p10.init)
          
          init.vals <- qlogis(c(psi.init, gamma.init, epsilon.init,
                                p11.init, p10.init, b.init))
        }
        
      }else{
        init.vals <- NULL
      }
      
      catch.error <- tryCatch(
        
        # Fit models for all params over each rep for this scenario
        m0 <- occMod(model = formulae,
                     data = one.pao,
                     type = model.type,
                     initvals = init.vals,
                     outfile = 'modname'),
        
        error = function(e) e
        
      ) # End tryCatch
      
      if (inherits(catch.error, "error")) {
        message(paste0(catch.error))
        
        # Predicted values are NAs because there is no prediction.
        occ.est[,i] <- col.est[,i] <-  ext.est[,i] <-
          p10.est[,i] <- p11.est[,i] <- b.est[,i] <- p.est[,i] <- NA
        
      }else{
        if (!(is.null(m0$warnings$conv)) | !(is.null(m0$warnings$VC))) {
          cw <- ifelse(!(is.null(m0$warnings$conv)),
                       cw <- as.character(m0$warnings$conv), cw <- NA)
          vw <- ifelse(!(is.null(m0$warnings$VC )), vw <- TRUE, vw <- FALSE)
          this.warn <- data.table(scenario = as.numeric(ex),
                                  rep = as.numeric(i),
                                  conv.warning = cw,
                                  VC.warning = vw)
          warns <- rbind(warns, this.warn)
        }
        
        # If no error, gather estimates ('predictions')
        
        # Gather predictions if using Miller model
        if (model.type == 'do.fp') {
          est.list <- lapply(m0$real, '[[', 1)
          occ.est[,i] <- est.list$psi
          col.est[,i] <- est.list$gamma
          ext.est[,i] <- est.list$epsilon
          p10.est[,i] <- est.list$p10
          p11.est[,i] <- est.list$p11
          b.est[,i] <- est.list$b
        }
        
        # Gather predictions if using regular dynamic model
        if (model.type == 'do.1') {
          occ.est[,i] <- m0$real$psi$est
          col.est[,i] <- m0$real$gamma$est
          ext.est[,i] <- m0$real$epsilon$est
          p.est[,i] <- m0$real$p$est
        }
        
        
      } # end if - else
      
    }
    
    # Compute the difference between estimated and actual
    diff.occ <- occ.est - occ.actuals
    diff.col <- col.est - col.actuals
    diff.ext <- ext.est - ext.actuals
    
    # Differences in detection parameters for Miller model
    if (model.type == 'do.fp') {
      diff.p10 <- p10.est - p10.actuals
      diff.p11 <- p11.est - p11.actuals
      diff.b <- b.est - b.actuals
    }
    
    # Differences in detection parameters for regular dynamic model
    if (model.type == 'do.1') {
      diff.p <- p.est - p11.actuals # p11 is equivalent to p in the regular model
    }
    
    # Gather state parameter bias info
    state.biases <- list(occ = diff.occ, col = diff.col, ext = diff.ext)
    b.state <- list()
    for (b in 1:length(state.biases)) {
      bias <- data.table::melt(data = state.biases[[b]],
                               id.vars = site,
                               measure.vars = paste0('rep.', 1:rep.scenario))
      bias$parameter <- names(state.biases)[[b]]
      b.state[[b]] <- bias
    }
    bias.state <- rbindlist(l = b.state)
    
    # Gather detection paramter bias info
    ifelse(model.type == 'do.fp',
           det.biases <- list(p11 = diff.p11, p10 = diff.p10, b = diff.b),
           det.biases <- list(p = diff.p))
    b.det <- list()
    for (b in 1:length(det.biases)) {
      dbias <- data.table::melt(data = det.biases[[b]],
                                id.vars = site,
                                measure.vars = paste0('rep.', 1:rep.scenario))
      dbias$parameter <- names(det.biases)[[b]]
      b.det[[b]] <- dbias
    }
    bias.detection <- rbindlist(l = b.det)
    
    # Remove "sites" since now we are doing an intercept model, so all site bias will be the same,
    # and there is no need to have all the duplicated rows
    bias.state[,Var1 := NULL]
    setkey(bias.state)
    unique.bias.state <- unique(bias.state)
    
    colnames(unique.bias.state) <- c('rep.scenario', 'bias', 'parameter')
    colnames(bias.detection) <- c('unit', 'rep.scenario', 'bias', 'parameter')
    unique.bias.scenario.state <- cbind(unique.bias.state, one.scen)
    bias.detection <- cbind(bias.detection, one.scen)
    
    # Save bias data
    bias <- list(bias.state = unique.bias.scenario.state,
                 bias.detection = bias.detection)
    saveRDS(bias, paste0(experiment.folder, 'bias_', model.name,
                         '_scenario_', scen, '.RDS'))
    
    # Collect and save any warnings
    saveRDS(object = warns, file = paste0(experiment.folder,'presence_warnings.RDS'))
  }
  
  cat('Bias calculations complete... gathering bias data.table.\n')
  
  # Collect biases into one big dt
  #   not saving the detection list biases for now -- too big and not sure if want them anyway
  bias.state.list <- bias.detection.list <-
    list(vector(mode = 'list', length = nrow(amdata@data$scenario.grid)))
  bias.files <- list.files(path = experiment.folder, pattern = 'bias')
  for (ex in seq(bias.files)) {
    which.file <- bias.files[ex]
    cat(ex, '...', which.file, '\n')
    
    # Files may not be in order; find the exact experiment name
    this.exp.number <- as.numeric(
      unlist(
        lapply(
          strsplit(which.file, split = '[_.]'), '[[', 4)
      )
    )
    
    # Assign contents to the amdatalist
    bias.scen <- readRDS(file = paste0(experiment.folder, which.file))
    bias.state.list[[ex]] <- bias.scen$bias.state
    if (save.detections) {
      
      # we are just going to take the AVERAGE bias for each parameter, for each rep.scenario (not all 6000 units ==> storage too large)
      bias.det.means <- bias.scen$bias.detection[,mean(bias), by = c('rep.scenario', 'parameter')]
      bias.det.means.mg <- merge(x = bias.det.means,
                                 y = bias.scen$bias.detection,
                                 by = c('rep.scenario', 'parameter'),
                                 all.x = TRUE)
      # Drop the unit column
      bias.det.means.mg[,c('unit', 'bias') := NULL]
      bias.det <- unique(bias.det.means.mg)
      colnames(bias.det)[which(colnames(bias.det) == 'V1')] <- 'bias' # replace with average bias
      bias.detection.list[[ex]] <- bias.det
    }
  }
  
  # Save state biases into one data.table (this should work fine unless thousands of sites)
  bias.state <- rbindlist(l = bias.state.list)
  saveRDS(object = bias.state, file = paste0(experiment.folder,'bias_state.RDS'))
  
  
  # Save AVERAGED detection biases into one table
  if (save.detections) {
    bias.detection <- rbindlist(l = bias.detection.list)
    saveRDS(object = bias.detection, file = paste0(experiment.folder,'bias_detection.RDS'))
  }
  
  # Return the state parameter bias data.table
  return(bias.state)
}






# simDetections ===============================================================
#' Simulate automated detection outcomes (TP:FN & TN:FP) for a focal species through time.
#' @name simDetections
#' @title Simulate automated detection outcomes for a focal species through time.
#' @param dynamics 'dynamics' object output from simDynamics()
#' @param vocalizations 'vocalizations' object output from simVocals()
#' @param false.alarm.rate False Positive rate (per minute). This value is as 'lambda' to draw poisson random False Positives in a given minute.
#' @param pr.TP.shape Shapes for pr.TP
#' @return List of data.tables with detections for each season.
#' @details tbd
#' @seealso Coming soon
#' @export
#' @examples
#'
#'\dontrun{
#' x <- NA
#'}

simDetections <- function(dynamics,
                          vocalizations,
                          false.alarm.rate,
                          pr.TP.shape
)
{
  
  # Extract helpful parameters for this function:
  n.seasons <- length(grep(x = names(dynamics), pattern = 'season'))
  n.sites <- length(unique(dynamics$site))
  n.days.in.season <- length(unique(vocalizations$season.day))
  daily.n.samples <- length(unique(vocalizations$recording))
  
  # Set prTp shape beta parameters:
  # (for journal submission, these are hardcoded in to reflect
  #  the beta paramters we used in the paper)
  if (pr.TP.shape == 'high_low') {
    tps <- c(4, 1)
    fps <- c(1, 4)
  }
  
  if (pr.TP.shape == 'med_med') {
    tps <- c(3, 3)
    fps <- c(3, 3)
  }
  
  # Loop through seasons to simulate automated detections:
  detections <- list()
  for (s in 1:n.seasons) {
    voc.season <- vocalizations[season == s]
    
    # Add a column for randomly generated False Positives:
    voc.season[, false.alarm := rpois(n = .N, lambda = false.alarm.rate)]
    
    # Subset voc.season to only keep rows with true.vocal and/or FP > 0
    voc.season.dets <- voc.season[true.vocal > 0 | false.alarm > 0]
    
    # Reshape to long
    dets.long <- data.table::melt(voc.season.dets,
                                  id.vars = c('site', 'season.day', 'recording'),
                                  measure.vars = c('true.vocal', 'false.alarm'))
    
    # Keep only rows where true vocalizations or FP detections occurred
    dets.long <- dets.long[value > 0]
    
    # Next, exlode the dets.long table out so that each "value" becomes its own observation
    dets.long <- dets.long[rep(seq(1, .N), value)]
    
    # Now that we have exploded it out, we can eliminate the value column
    dets.long[, value := NULL]
    
    # Next, we assign the pr.TP for each detected event
    #   true.vocals have their pr.TP drawn from the pr.TP.shapes beta distribution
    #   false.alarms have their pr.TP drawn from the pr.FP.shapes beta distribution
    dets.long[variable == 'true.vocal',
              pr.TP := rbeta(n = .N, shape1 = tps[1], shape2 = tps[2])]
    dets.long[variable == 'false.alarm',
              pr.TP := rbeta(n = .N, shape1 = fps[1], shape2 = fps[2])]
    colnames(dets.long)[colnames(dets.long) == 'variable'] <- 'vocal.type'
    setkeyv(x = dets.long, cols = c('site', 'season.day', 'recording'))
    
    # Save detections from this season
    dets.long[ ,season := s]
    
    # Add presence/absence status in order to be able to calculate
    #   actual singing rates (given presence) later on
    pres <- dynamics[, c('site', paste0('season.',s)), with = FALSE]
    names(pres)[2] <- 'present'
    merge.dets <- merge(x = dets.long, y = pres, by = 'site', all.x = TRUE)
    detections[[s]] <- merge.dets
    names(detections)[[s]] <- paste0('season.', s)
    
  } # end for season
  
  detections.dt <- rbindlist(l = detections)
  
}





# simDynamics =================================================================
#' Simulate dynamics of underlying/latent/unobservable occurrence, colonization, extinction, and vocalization patterns through time for a focal species.
#' @name simDynamics
#' @title Simulate underlying patterns of occurrence, colonization, and extinction patterns for focal species.
#' @param occupancy.prob Probability of occupancy
#' @param colonization.prob Probility of colonization, given absence
#' @param extinction.prob Probability of extinction, given presence
#' @param n.sites Number of sites to include in the simulation.
#' @param n.seasons Number of seasons to include in the simulation.
#' @param sampling.type Type of site selection sampling to perform. 'random' provides a random sample. 'stratified' provides a stratified sample, with options to pass stratify.by and n.groups arguments to the stratifySiteSample function. NOTE: for journal submission, this argument is a default of 'random' and will not work otherwise.
#' @param ... Additional arguments fed into stratifySiteSample()
#' @return A data table of occurrence dynamics
#' @details tbd...
#' @seealso Coming soon
#' @export
#' @examples
#'\dontrun{
#' x <- NA
#'}
#'

simDynamics <- function(occupancy.prob,
                        colonization.prob,
                        extinction.prob,
                        n.sites,
                        n.seasons,
                        sampling.type = 'random',
                        ...
)
{
  # Set up a data.table to collect occurrence dynamics probabilities
  seasons <- matrix(0, nrow = n.sites, ncol = n.seasons)
  dyn <- data.table(cbind(1:n.sites,
                          occupancy.prob,
                          colonization.prob,
                          extinction.prob,
                          seasons))
  colnames(dyn) <- c('site', 'pr.occ', 'pr.col', 'pr.ext',
                     paste0('season.', 1:n.seasons))
  dyn$season.1 <- rbinom(n = n.sites, size = 1, prob = dyn[['pr.occ']])
  season.columns <- grep(x = names(dyn), pattern = 'season')
  
  # Simulate occurrence dynamics for each season
  for (i in 2:n.seasons) { # start with the second season
    # calculate occurrence based on ifs.
    last.season.column <- season.columns[i - 1]
    last.season <- dyn[,..last.season.column]
    this.season.column <- season.columns[i]
    pr.ext <- dyn$pr.ext[last.season == 1]
    pr.col <- dyn$pr.col[last.season == 0]
    dyn[which(last.season == 1), (this.season.column) :=
          as.numeric(rbinom(n = length(pr.ext), size = 1, prob = 1 - pr.ext))]
    dyn[which(last.season == 0), (this.season.column) :=
          as.numeric(rbinom(n = length(pr.col), size = 1, prob = pr.col))]
  }
  setkey(dyn, site)
  
  # Return a data.table of simulated occurrence dynamics
  dyn
}

# simOccupancyExperiment ======================================================
#' Conduct a simulated multi-season automated acoustic monitoring experiment.
#' @name simOccupancyExperiment
#' @title Conduct a simulated multi-season automated acoustic monitoring experiment.
#' @param scenario.grid Grid of scenarios to run produced by simScenarioGrid()
#' @param experiment.folder Path to store RDS results of each scenario. This function sends each scenario to a folder separately for each scenario, ensuring that each scenario has been stored in case an edge case has found and the simulation fails for some reason. This is a temporary convenience for Cathleen in dev phase, and may change in the final version of AMMonitor software. 
#' @param voc.time.units Time units vocalization probabilities or rates were modeled on. 'minutely' if probabilities or rates should be applied to each minute. 'hourly' if probabilities or rates should be applied to each hour.
#' @param rep.scenario How times to repeat a single scenario
#' @param ... Additional arguments to simDynamics, stratifySiteSample.
#' @return Saves amdata objects to the specified experiment.folder. 
#' @details tbd.
#' @seealso Coming soon
#' @export
#' @examples
#'\dontrun{
#'x <- NA
#'}
#'

simOccupancyExperiment <- function(
  scenario.grid,
  experiment.folder,
  voc.time.units,
  rep.scenario,
  ...
)
  
{
  
  # Add a forward slash to experiment folder path if missing
  if (grepl("\\/$", experiment.folder) == FALSE) {
    experiment.folder <- paste0(experiment.folder, '/')
  }
  
  # Sort the parameter scenario grid
  setkey(scenario.grid, n.sites, daily.n.samples, n.seasons,
         occ, col, ext, voc.rate, any.true, n.days)
  
  cat('This experiment contains', nrow(scenario.grid), 'scenarios. \n')
  cat('Each scenario will be run', rep.scenario, 'times. \n')
  
  #  Run experiment
  for (ex in 1:nrow(scenario.grid)) {
    
    scen <- scenario.grid$scenario[ex]
    cat('This is scenario', scen,'\n')
    
    # Set up parameters for experiment 'ex'
    occ.prob <- scenario.grid[scenario == scen, occ]
    col.prob <- scenario.grid[scenario == scen, col]
    ext.prob <- scenario.grid[scenario == scen, ext]
    voc.rate <- scenario.grid[scenario == scen, voc.rate]
    n.seasons <- scenario.grid[scenario == scen, n.seasons]
    n.days.in.season <- scenario.grid[scenario == scen, n.days.in.season]
    daily.n.samples <- scenario.grid[scenario == scen, daily.n.samples]
    pr.TP.vals <- scenario.grid[scenario == scen, pr.TP]
    n.sites <- scenario.grid[scenario == scen, n.sites]
    pct.conf <- scenario.grid[scenario == scen, pct.confirm]
    false.alarm.rate <- scenario.grid[scenario == scen, false.alarm.rate]
    n.days <- scenario.grid[scenario == scen, n.days]
    any.true <- scenario.grid[scenario == scen, any.true]
    cull.FP <- scenario.grid[scenario == scen, cull.FP]
    
    # Run iterations for this experiment:
    dyn.list <- enc.list <- list()
    sound.summary <- matrix(data = 0, nrow = rep.scenario, ncol = 5)
    colnames(sound.summary) <- c('lambda.c.presence', 'lambda.c.overall',
                                 'lambda.f', 'true.pr.TP', 'false.pr.TP')
    for (i in 1:rep.scenario) {
      
      cat(i, '\n')
      
      # Simulate dynamics:
      dynamics <- simDynamics(occupancy.prob = occ.prob,
                              colonization.prob = col.prob,
                              extinction.prob = ext.prob,
                              n.sites = n.sites,
                              n.seasons = n.seasons,
                              ...)
      
      # Simulate vocalizations:
      vocs <- simVocals(vocalization.rate = voc.rate,
                        voc.time.units = voc.time.units,
                        dynamics = dynamics,
                        n.days.in.season = n.days.in.season,
                        daily.n.samples = daily.n.samples)
      
      # Simulate detections
      dets <- simDetections(dynamics = dynamics,
                            vocalizations = vocs,
                            pr.TP.shape = pr.TP.vals,
                            false.alarm.rate = false.alarm.rate)
      
      # Aggregate the surveys
      enc.hists <- simAggregate(dynamics = dynamics,
                                detections = dets,
                                vocalizations = vocs,
                                pct.confirm = pct.conf,
                                n.days.in.survey = n.days,
                                any.true.threshold = any.true,
                                cull.FP = cull.FP)
      
      # Save summary results for detections (this occurs for reporting results in journal submission)
      
      # Number of sites with presence (summed across both seasons)
      n.sites.present <- sum(dets[present == 1, length(unique(site)), by = season]$V1)
      
      # Calculate singing rate result (given presence)
      lambda.c.presence <- dets[present == 1 & vocal.type == 'true.vocal', .N] /
        (n.sites.present*n.days.in.season*daily.n.samples)
      
      # ^we don't need to multiply by n.seasons here because we are summing the sites already
      lambda.c.overall <- dets[vocal.type == 'true.vocal', .N] /
        (n.sites*n.days.in.season*daily.n.samples)
      
      # Calculate false alarm rate results
      lambda.f <- dets[vocal.type == 'false.alarm', .N] /
        (n.sites*n.days.in.season*n.seasons*daily.n.samples)
      
      # mean prTPs (will be same regardless of presence)
      mean.pr.TP <- dets[,mean(pr.TP), by = vocal.type]
      sound.summary[i, ] <- c(lambda.c.presence, lambda.c.overall,lambda.f,
                              mean.pr.TP[vocal.type == 'true.vocal', V1],
                              mean.pr.TP[vocal.type == 'false.alarm', V1])
      
      dyn.list[[i]] <- dynamics
      enc.list[[i]] <- enc.hists
      
    } # end 'i' rep.scenario loop for each scenario
    
    # Save each scenario separately as an RDS first to avoid loop breakage
    single.scen <- list(scenario.id = scen,
                        sound.summary = sound.summary,
                        dynamics = dyn.list,
                        encounter.histories = enc.list)
    saveRDS(object = single.scen,
            file = paste0(experiment.folder, 'scenario_', scen,'.RDS'))
    
  } # end 'ex' parallel loop
  
  message('Experiment complete.')

  # Collect all experiment results into amdata
  simScenarioCollect(
    experiment.folder = experiment.folder,
    scenario.grid = scenario.grid
  )
  
} # end function

# simPlotBias =================================================================
#' Make dotplots of bias for manuscript figures.
#' @name simPlotBias
#' @title Make dotplots of bias for manuscript figures.
#' @param experiment.folder Path to store RDS results of each experiment. 
#' @param false.alarm.rate.paper A SINGLE VALUE of the false alarm rate to use for these plots. (e.g. 0.8)
#' @param pct.confirms.paper Vector of pct.confirms to use in the paper (e.g. c(0.025, 0.5)). ONLY PUT TWO IN AT A TIME. otherwise the figures become overwhelming.
#' @param any.true.thresholds.paper Vector of any.true thresholds to use in the paper (e.g. c(0.8, 0.95)) 
#' @param n.days.paper Vector (of no more than 2) n.days to use in the paper (e.g. c(1,3)). Currently, this function expects 1 and/or 3 day aggregation only.
#' @param plot.state.parameters Logical flag for which sets of params to plot. If TRUE (default) plots bias for the state parameters. If FALSE, plots bias for the detection parameters.
#' @param remove.warnings Remove experiment replicates that had warnings. Options: 'none' removes no replicates, 'vc' only removes replicates that encountered a warning in the variance-covariance matrix, 'convergence' only removes replicates that encountered convergence warnings, and 'both' removes any replicate that had either a vc or a convergence warning. Default = 'convergence'.
#' @param zoom Logical flag for whether to zoom
#' @param legend Logical flag for whether the plot is only being used to generate the legend. Default false.
#' @param scale.y.continuous ggplot scale_y_breaks object to set the y-axis limits and the y-axis breaks. Default = scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.25), limits = c(-1, 1))
#' @return awesome ggplots
#' @details This is a function created to make plotting easier. I provide no guarantee it will work outside of the argument inputs used for the journal submission experiment. 
#' @seealso Coming soon
#' @export
#' @examples
#'
#' #'/dontrun{
#'x <- NA
#'}
#

simPlotBias <- function(experiment.folder,
                        false.alarm.rate.paper,
                        any.true.thresholds.paper,
                        n.days.paper,
                        pct.confirms.paper,
                        plot.state.parameters = TRUE,
                        remove.warnings = 'convergence',  #c('both', 'convergence', 'vc', 'none'),
                        zoom = FALSE,
                        legend = FALSE,
                        scale.y.continuous = scale_y_continuous(
                          breaks = seq(from = -1, to = 1, by = 0.25),
                          limits = c(-1, 1)))
{

  # Add a forward slash to experiment folder path if missing
  if (grepl("///$", experiment.folder) == FALSE) {
    experiment.folder <- paste0(experiment.folder, '/')
  }
  
  if (legend) zoom <- TRUE
  
  # Read in necessary files:
  ifelse(plot.state.parameters == TRUE,
         bias <- readRDS(paste0(experiment.folder,'bias_state.RDS')),
         bias <- readRDS(paste0(experiment.folder,'bias_detection.RDS')))
  warns <- readRDS(paste0(experiment.folder,'presence_warnings.RDS'))
  experiment <- readRDS(paste0(experiment.folder,'amdata_object.RDS'))
  s.grid <- experiment@data$scenario.grid
  s.grid.paper <- s.grid[pct.confirm %in% pct.confirms.paper &
                           n.days %in% n.days.paper &
                           false.alarm.rate == false.alarm.rate.paper]
  warns <- warns[scenario %in% s.grid.paper$scenario]
  
  # Capture convergence rate
  rep.scenario <- length(experiment@data$encounter.histories[[1]])
  convergence.rate <- data.table(1 - table(warns$scenario)/rep.scenario)
  colnames(convergence.rate) <- c('scenario', 'convergence.rate')
  convergence.rate$scenario <- as.numeric(convergence.rate$scenario)
  setkey(convergence.rate)
  convergence.rate <- unique(convergence.rate)
  bias <- merge(x = convergence.rate, y = bias, by = 'scenario', all.x = TRUE)
  bias[is.na(convergence.rate), convergence.rate := 1] # if NA, that means there were no warnings, 100% converge
  
  # Add data about rep warnings so that unconverged reps can be eliminated from viz
  warns[,rep.scenario := paste0('rep.', rep)]
  
  if (remove.warnings == 'none') {
    # Treat all warned scenario replicates as if they have converged
    warns[,converged := TRUE]
  }
  
  if (remove.warnings == 'vc') {
    # Treat scenario replicates with a VC warning as unconverged
    warns[VC.warning == TRUE, converged := FALSE]
    warns[VC.warning == FALSE, converged := TRUE]
  }
  
  if (remove.warnings == 'convergence') {
    # Treat scenario replicates with a conv.warning as unconverged
    warns[is.na(conv.warning), converged := TRUE]
    warns[!(is.na(conv.warning)), converged := FALSE]
  }
  
  if (remove.warnings == 'both') {
    # Treat all warned scenario replicates as if they have not converged
    warns[,converged := FALSE]
  }
  
  # Merge warnings with bias table for plotting
  b <- merge(x = bias,
             y = warns[,c('scenario', 'rep.scenario', 'converged')],
             by = c('scenario', 'rep.scenario'),
             all.x = TRUE)
  b[is.na(converged), converged := TRUE] # reps w/no warns have converged
  
  # Prep data.table for plotting
  b[occ == 0.6 & col == 0.25 & ext == 0.25, state := 'HH']
  b[occ == 0.2 & col == 0.25 & ext == 0.25, state := 'LH']
  b[occ == 0.6 & col == 0.05 & ext == 0.05, state := 'HL']
  b[occ == 0.2 & col == 0.05 & ext == 0.05, state := 'LL']
  b[voc.rate == 100, state.voc := paste0(state,'H')]
  b[voc.rate == 20, state.voc := paste0(state,'L')]
  b[voc.rate == 100, vr := 'High Call Rate']
  b[voc.rate == 20, vr := 'Low Call Rate']
  b[col == 0.25, turnover := 'high']
  b[col == 0.05, turnover := 'low']
  b[pr.TP == 'high_low', classifier := 'Good Classifier']
  b[pr.TP == 'med_med', classifier := 'Bad Classifier']
  
  # Only keep converged replicates:
  bias <- b[pct.confirm %in% pct.confirms.paper &
              n.days %in% n.days.paper &
              false.alarm.rate == false.alarm.rate.paper
            & converged == TRUE]
  
  # Compute the number of figures per panel
  n.figs.per.panel <- length(pct.confirms.paper)*length(unique(n.days.paper))
  
  # set vector of pct.confirms to loop through -- sorted
  pc <- sort(rep(pct.confirms.paper, length(pct.confirms.paper)))
  
  # set vector of ndays aggregation to loop through -- unsorted
  nd <- rep(unique(bias$n.days), length(unique(bias$n.days)))
  
  if (zoom) {
    pc <- unique(pc)
    nd <- unique(nd)
  }
  
  # total number of figures
  n.figs.total <- n.figs.per.panel*length(any.true.thresholds.paper)
  
  # set titles: (hardcoded these in for my own convenience - cb)
  if (zoom) {
    fig.titles <- rep('1-day Aggregation (30 Surveys)', 2)
    zoom.tag <- 'zoom' # create a tag for the figure filepath
  }
  
  if (zoom == FALSE) {
    fig.titles <- c('1-day Aggregation (30 Surveys)', '3-day Aggregation (10 Surveys)')
    fig.letters <- letters[1:n.figs.total]
    fig.titles <- rep(fig.titles, n.figs.per.panel/2)
    fig.titles <- paste0(fig.letters, '. ', fig.titles)
    zoom.tag <- 'full' # create a tag for the figure filepath
  }
  
  # Set general themes and niceties:
  ifelse(legend == TRUE, leg.pos <- 'left', leg.pos <- 'none')
  mytheme <- theme(axis.ticks.y = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(color = 'black'),
                   text = element_text(family = 'Times'),
                   plot.title = element_text(hjust = 0.5, size = 10),
                   legend.text = element_text(size = 10),
                   legend.position = leg.pos,
                   legend.title = element_blank(),
                   strip.background = element_blank(),
                   axis.title.x = element_text(size = 10))
  windowsFonts(Times = windowsFont('Times New Roman'))
  dodge <- position_dodge(0.75)
  bias$state <- factor(bias$state) # levels: HH 1, HL 2, LH 3, LL 4
  
  ifelse(plot.state.parameters,
         bias$parameter <- factor(bias$parameter, levels = c('occ', 'col', 'ext')),
         bias$parameter <- factor(bias$parameter, levels = c('p11', 'p10', 'b')))
  
  # Create a tag for the figure filepath
  ifelse(plot.state.parameters, param.tag <- 'state', param.tag <- 'detection')
  
  # Set appropriate labels
  if (plot.state.parameters) {
    scale_x_discrete_values <-
      scale_x_discrete(labels = c('occ' = expression(psi),
                                  'col' = expression(gamma),
                                  'ext' = expression(epsilon)))
    xlabel <- xlab('State Parameter')
  }else{
    scale_x_discrete_values <- scale_x_discrete(labels = c('p11', 'p10', 'b'))
    xlabel <- xlab('Detection Parameter')
  }
  
  # Loop through any.true.thresholds to create mutiple windows:
  for (gr in 1:length(any.true.thresholds.paper)) {
    
    any.true.threshold <- any.true.thresholds.paper[[gr]]
    
    gg <- list()
    
    for (i in 1:n.figs.per.panel) {
      
      if (zoom) {
        these.titles <- fig.titles
        ndi <- nd
      }else{
        if (length(any.true.thresholds.paper) == 1) {
          these.titles <- fig.titles
        }else{
          ifelse(gr == 1,
                 these.titles <- fig.titles[(n.figs.per.panel + 1):n.figs.total],
                 these.titles <- fig.titles[1:n.figs.per.panel])
        }
        ndi <- nd[i]
        
        # # check if there are more than 1 pct.confirm
        # ifelse(length(pc) == 1, pci <- pc, pci <- pc[i])
      }
      
      # check if there are more than 1 pct.confirm
      ifelse(length(pc) == 1, pci <- pc, pci <- pc[i])
      
      dt <- bias[false.alarm.rate == false.alarm.rate.paper &
                   any.true == any.true.threshold &
                   n.days == ndi &
                   converged == TRUE &
                   pct.confirm == pci]
      
      # summarize means and sds for each level:
      mns <- dt[, .(mean(bias), sd(bias)),
                by = c('state', 'classifier', 'vr', 'parameter')]
      colnames(mns)[5:6] <- c('mean.bias', 'sd.bias')
      
      # create the summary table, keeping only columns needed for plotting
      dt2 <- merge(x = dt[,c('state', 'classifier', 'vr', 'parameter')],
                   y = mns, by = c('state', 'classifier', 'vr', 'parameter'))
      setkey(dt2)
      dt <- unique(dt2)
      dt[vr == 'High Call Rate' & classifier == 'Bad Classifier',
         classif.vr := 'Bad Classifier + High Call Rate']
      dt[vr == 'Low Call Rate' & classifier == 'Bad Classifier',
         classif.vr := 'Bad Classifier + Low Call Rate']
      dt[vr == 'High Call Rate' & classifier == 'Good Classifier',
         classif.vr := 'Good Classifier + High Call Rate']
      dt[vr == 'Low Call Rate' & classifier == 'Good Classifier',
         classif.vr := 'Good Classifier + Low Call Rate']
      
      # Plot
      gg[[i]] <- ggplot(data = dt, aes(x = parameter,
                                       y = mean.bias,
                                       color = classif.vr,
                                       shape = classif.vr)) +
        geom_point(position = dodge) +
        geom_errorbar(aes(ymin = mean.bias - sd.bias, ymax = mean.bias + sd.bias),
                      position = dodge,
                      width = 0.05,
                      linetype = 'dotted') +
        facet_grid( ~ state) +
        scale_color_manual(values = c('gray50', 'gray50', 'black', 'black')) +
        scale_shape_manual(values = c(16, 1, 16, 1)) +
        scale.y.continuous +
        scale_x_discrete_values +
        xlabel +
        ggtitle(these.titles[i]) +
        ylab('Bias') +
        theme_classic() +
        mytheme
    } # end i n.figs.per.panel
    
    # Height of figure needs to automatically scale with number of pct.confirms
    ht <- length(unique(pc))*3
    
    if (zoom) {
      x11(height = ht, width = 5.5)
      
      first.title <- textGrob(paste0('a. ', unique(pc)[1]*100,
                                     ' % Confirmation, Survey-level Detection Threshold = ',
                                     any.true.threshold,
                                     ', N Days = ', unique(nd)),
                              gp = gpar(fontfamily = 'Times'))
      second.title <- textGrob(paste0('b. ', unique(pc)[2]*100,
                                      ' % Confirmation, Survey-level Detection Threshold = ',
                                      any.true.threshold,
                                      ', N Days = ', unique(nd)),
                               gp = gpar(fontfamily = 'Times'))
      first.g <- grid.arrange(gg[[1]], ncol = 1, top = first.title)
      second.g <- grid.arrange(gg[[2]], ncol = 1, top = second.title)
      pl <- grid.arrange(first.g, second.g)
      
    }else{
      x11(height = ht, width = 7)
      
      super.title <- textGrob(paste0(unique(pc)[1]*100,
                                     ' % Confirmation, Survey-level Detection Threshold = ',
                                     any.true.threshold),
                              gp = gpar(fontfamily = 'Times'))
      first.g <- grid.arrange(gg[[1]], gg[[2]], ncol = 2, top = super.title)
      
      if (length(gg) > 2) { # only will plot 2-4 at a time. Only put in 2 pct.confirms at a time.
        super.title.2 <- textGrob(paste0(unique(pc)[2]*100,
                                         ' % Confirmation, Survey-level Detection Threshold = ',
                                         any.true.threshold),
                                  gp = gpar(fontfamily = 'Times'))
        second.g <- grid.arrange(gg[[3]], gg[[4]], ncol = 2, top = super.title.2)
        pl <- grid.arrange(first.g, second.g)
        
      }
    } # end if not zoom
  } #end for any.true.thresholds grid (gr)
} # end loopy function


# simScenarioCollect ==========================================================
#' Collect individual scenario RDS files into a single amdata object
#' @name simScenarioCollect
#' @title Collect individual scenario RDS files into a single amdata object
#' @param scenario.grid Grid of scenarios to run produced by simScenarioGrid()
#' @param experiment.folder Path to store RDS results of each scenario 
#' @return Returns a single Amdata objects into a folder.
#' @details tbd...
#' @seealso Coming soon
#' @export
#' @examples
#'
#' \dontrun{
#' x <- NA
#' }
#'

simScenarioCollect <- function(experiment.folder,
                               scenario.grid){
  
  # Add a forward slash to experiment folder path if missing
  if (grepl("\\/$", experiment.folder) == FALSE) {
    experiment.folder <- paste0(experiment.folder, '/')
  }
  
  # Read in all the RDS experiment files:
  experiment.files <- list.files(path = experiment.folder, pattern = '^scenario')
  
  # Preallocate the amdatalist
  amdatalist <- list(dynamics = vector(mode = 'list',
                                       length = length(experiment.files)), # for psi, e, c
                     encounter.histories = vector(mode = 'list',
                                                  length = length(experiment.files)))
  sound.summaries <- vector(mode = 'list',
                            length = length(experiment.files))
  
  cat('Collecting results...\n')
  
  for (ex in seq(experiment.files)) {
    
    which.file <- experiment.files[ex]
    
    cat(ex, '...', which.file, '\n')
    
    # Files may not be in order; find the exact experiment name
    this.scen.number <- as.numeric(
      unlist(
        lapply(
          strsplit(which.file, split = '[_.]'), '[[', 2)
      )
    )
    # Assign contents to the amdatalist
    scen <- readRDS(file = paste0(experiment.folder, which.file))
    sound.sum <- data.table(scen$sound.summary)
    sound.sum[,scenario := this.scen.number]
    amdatalist$dynamics[[this.scen.number]] <- scen$dynamics
    amdatalist$encounter.histories[[this.scen.number]] <- scen$encounter.histories
    sound.summaries[[this.scen.number]] <- sound.sum
  }
  
  # Create detection summaries (calculate p11, p10, b)
  rep.scenario <- length(amdatalist$dynamics[[1]])
  n.seasons <- length(grep(colnames(amdatalist$dynamics[[1]][[1]]), pattern = 'season'))
  n.sites <- nrow(amdatalist$dynamics[[1]][[1]])
  survey.fp.list <- list()
  
  for (i in 1:length(amdatalist$encounter.histories)) {
    
    # Create a matrix to characterize the following at the survey-level:
    #   - uncertain 1s that were incorrect (false positives)
    #   - uncertain 1s that were correct
    #   - certain trues
    survey.fp.scenario <- matrix(0, ncol = 3, nrow = rep.scenario)
    survey.fp.scenario <- data.table(survey.fp.scenario)
    survey.fp.scenario[,scenario := i]
    names(survey.fp.scenario)[1:3] <- c('p10', 'p11', 'b')
    
    # For each scenario replicate
    for (j in 1:rep.scenario) {
      
      # Extract dynamics and ecounter histories 
      dyn <- amdatalist$dynamics[[i]][[j]]
      eh <- amdatalist$encounter.histories[[i]][[j]]
      
      survey.seasons <- matrix(0, ncol = 3, nrow = 2)
      names(survey.seasons) <- c('p10', 'p11', 'b')
      
      for (s in 1:n.seasons) {
        dyn.s <- dyn[, paste0('season.',s), with = FALSE]
        eh.s <- eh[, grep(x = colnames(eh), pattern = paste0(s,'-'))]
        
        occupied <- eh.s[which(dyn.s == 1), ]   # eh for occupied sites
        unoccupied <- eh.s[which(dyn.s == 0), ] # eh for unoccupied sites
        
        total.fp.1s <- sum(unoccupied)
        total.utp.1s <- sum(occupied == 1)
        total.ctp.2s <- sum(occupied == 2)
        
        ratio.fp <- total.fp.1s/(ncol(unoccupied)*nrow(unoccupied))
        ratio.utp <- total.utp.1s/(ncol(occupied)*nrow(occupied))
        ratio.ctp <- total.ctp.2s/(ncol(occupied)*nrow(occupied))
        
        survey.seasons[s, 1:3] <- c(ratio.fp, ratio.utp, ratio.ctp)
      } # end calc season
      vals <- colMeans(survey.seasons)
      survey.fp.scenario[j, p10 := vals[1]]
      survey.fp.scenario[j, p11 := vals[2]]
      survey.fp.scenario[j, b := vals[3]]
    } # end rep scenario
    survey.fp.list[[i]] <- survey.fp.scenario
  } # end this scenario
  
  detection.summary.dt <- rbindlist(l = survey.fp.list)
  amdatalist$detection.summary <- detection.summary.dt
  
  # Store the sound summaries for easy journal results reporting
  sound.summaries.dt <- rbindlist(l = sound.summaries)
  amdatalist$sound.summary <- sound.summaries.dt
  
  # Store the scenario grid containing all param combos:
  amdatalist$scenario.grid <- scenario.grid
  
  # Save as amdata object:
  amdata <- amData(data = amdatalist, comment = 'experiment')
  saveRDS(amdata, paste0(experiment.folder, 'amdata_object.RDS'))
  cat('Data saved as',  paste0(experiment.folder, 'amdata_object.RDS'))
}

# simScenarioGrid =============================================================
#' Create a grid of scenarios to use in simOccupancyExperiment()
#' @name simScenarioGrid
#' @title Create a grid of scenarios to use in simOccupancyExperiment()
#' @param occupancy.probs Probability of occupancy
#' @param colonization.probs Probility of colonization, given absence
#' @param extinction.probs Probability of extinction, given presence
#' @param n.seasons Number of seasons
#' @param n.days.in.season Number of days in each a season (assumed same n.days.in.season for each season in a single experiment.)
#' @param n.sites Numeric vector of the number of sites to simulate.
#' @param daily.n.samples How many one-minute recording samples to take per day. Assume they occur at ideal sampling times.
#' @param vocalization.rates lambda rate(s) of true target signal to use in rpois
#' @param false.alarm.rates lambda rate(s) of false alarm generative sources to use in rpois
#' @param pr.TP.shapes vector of distributions. E.g. c('high_low', 'med_med'). This would indicate that in the first scenario, prTP is high for true positives, and low for false positives (a good classifier). In the second scenario, prTP is medium for true positives, medium for false positives (a  bad classiifer). For simplicity, the actual beta values get constructed from hardcoded numbers later on in the simulation. Not generalizable right now because I don't need it to be.
#' @param cull.FP Logical vector of whether remove detections with target signal probabilities below 0.5. If TRUE, all detections with prTP < 0.5 will be culled from the dataset before aggregation.
#' @param any.true.thresholds Numeric vector of probability aggregation thresholds above which a True Positive is concluded.
#' @param n.days.in.survey Numeric vector of the number of days over which to aggregate TP probabiltiies.
#' @param pct.confirm Numeric vector of percentage of surveys to confirm (e.g. c(0.05, 1, 2))
#' @return Amdata objects into a folder...?
#' @details tbd...
#' @seealso Coming soon
#' @export
#' @examples
#'
#' \dontrun{
#' x <- NA
#' }
#'



simScenarioGrid <- function(
  
  # Latent models:
  occupancy.probs,
  colonization.probs,
  extinction.probs,
  vocalization.rates,
  
  # Season parameters:
  n.seasons,
  n.days.in.season,
  
  # Additional Factorial Dials:
  
  # ==> Dynamics Dials:
  n.sites,
  
  # ==> Miller et al. 2013 Dials
  pct.confirm,
  
  # ==> Detection Dials:
  pr.TP.shapes,
  daily.n.samples,
  false.alarm.rates,
  
  # ==> Aggregation Dials:
  any.true.thresholds,
  n.days.in.survey,
  cull.FP
  
){
  
  s.grid <- data.table(expand.grid(occ = occupancy.probs,
                                   col = colonization.probs,
                                   ext = extinction.probs,
                                   voc.rate = vocalization.rates,
                                   n.seasons = n.seasons,
                                   n.days.in.season = n.days.in.season,
                                   n.sites = n.sites,
                                   pct.confirm = pct.confirm,
                                   pr.TP = pr.TP.shapes,
                                   daily.n.samples = daily.n.samples,
                                   false.alarm.rate = false.alarm.rates,
                                   n.days = n.days.in.survey,
                                   any.true = any.true.thresholds,
                                   cull.FP = cull.FP,
                                   stringsAsFactors = FALSE))
  
  setkey(s.grid, n.sites, daily.n.samples, n.seasons,
         occ, col, ext, voc.rate, any.true, n.days, pct.confirm, cull.FP)
  
  s.grid[,scenario := 1:.N]
  
  # Only keep cases that fit the Miller 2015 specifcation of high occurrence, low turnover etc.
  #  (this aspect of the function is strictly for the journal paper)
  c1 <- s.grid[occ == 0.2 & col == 0.05 & ext == 0.05, scenario]   # low occurrence, low turnover
  c2 <- s.grid[occ == 0.2 & col == 0.25 & ext == 0.25, scenario]   # low occurrence, high turnover
  c3 <- s.grid[occ == 0.6 & col == 0.05 & ext == 0.05, scenario]   # high occurrence, low turnover
  c4 <- s.grid[occ == 0.6 & col == 0.25 & ext == 0.25, scenario]   # high occurrence, high turnover
  
  narrow.s.grid <- s.grid[scenario %in% c(c1, c2, c3, c4)]
  
  narrow.s.grid[,scenario := 1:.N] # correct the exp.names
  
  cat('This experiment contains', nrow(narrow.s.grid), 'scenarios. \n')
  
  narrow.s.grid
}

# simVocals ===================================================================
#' Simulate underlying vocalization patterns through time for a focal species.
#' @name simVocals
#' @title Simulate underlying vocalization patterns for a focal species.
#' @param vocalization.rate lambda rate of true target signal to use in rpois
#' @param voc.time.units Time units vocalization probabilities or rates were modeled on. 'minutely' if probabilities or rates should be applied to each minute. 'hourly' if probabilities or rates should be applied to each hour.
#' @param dynamics 'Dynamics' object output from simDynamics()
#' @param n.days.in.season Number of days in each a season (assumed same n.days.in.season for each season in a single experiment.)
#' @param daily.n.samples How many one-minute recording samples to take per day. Assume they occur at ideal sampling times.
#' @return A data table of awesome
#' @details tbd...
#' @seealso Coming soon
#' @export
#' @examples Coming soon
#'
#'\dontrun{
#' x <- 'example coming soon'
#'}

simVocals <- function(vocalization.rate,
                      voc.time.units,
                      dynamics,
                      n.days.in.season,
                      daily.n.samples
)
{
  
  if (missing(voc.time.units)) {
    stop('Please use the \'voc.time.units\' argument to clarify whether your modeled vocalization probabilities or rates should be applied on an hourly or minutely basis.')
  }
  
  # Extract necessary param info for the function:
  n.seasons <- length(grep(x = names(dynamics), pattern = 'season'))
  n.sites <- length(unique(dynamics$site))
  seas.cols.dyn <- grep(pattern = 'season', x = names(dynamics))
  nrows <- daily.n.samples*n.days.in.season*n.sites
  ifelse(voc.time.units == 'hourly',
         voc.rate <- vocalization.rate/60,
         voc.rate <- vocalization.rate)
  
  # Simulate underlying vocalization patterns for each season
  voc.seasons <- list()
  for (s in 1:n.seasons) {
    
    # Perform poisson draw for vocalizations for each site and sampling occasion:
    vocals <- data.table(site = sort(rep(1:n.sites, n.days.in.season*daily.n.samples)),
                         season.day = rep(sort(rep(1:n.days.in.season,
                                                   daily.n.samples)), n.sites),
                         recording = rep(1:daily.n.samples, n.sites*n.days.in.season),
                         true.vocal = rpois(n = nrows, lambda = voc.rate))
    
    # If species is absent from site, correct true vocalizations to 0:
    for (which.site in 1:n.sites) {
      seas.dyn <- seas.cols.dyn[s]
      status <- dynamics[site == which.site, ..seas.dyn]
      if (status == 0) {
        vocals[site == which.site, true.vocal := 0]
      }
    }# end presence 0 correction to vocalizations
    
    # Add presence/absence status in order to be able to calculate
    # #   actual singing rates (given presence) later on
    pres <- dynamics[, c('site', paste0('season.',s)), with = FALSE]
    names(pres)[2] <- 'present'
    merge.vocs <- merge(x = vocals, y = pres, by = 'site', all.x = TRUE)
    merge.vocs[,season := s]
    voc.seasons[[s]] <- merge.vocs
    names(voc.seasons)[[s]] <- paste0('season.', s)
    
  } # end seasons loop
  
  # Return data table of simulated vocalizations
  voc.seasons.dt <- rbindlist(l = voc.seasons)
  
}
