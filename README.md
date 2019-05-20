Code to reproduce analysis from "Dynamic wildlife occupancy models using automated acoustic monitoring data" (Balantic &amp; Donovan 2019, Ecological Applications, DOI: 10.1002/eap.1854)

This Github repository exists in service of open science, to reproduce the analysis from the following paper:

Balantic, C. M., & Donovan, T. M. (2019).
Dynamic wildlife occupancy models using automated acoustic monitoring data. Ecological Applications. https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/eap.1854

Repository code contained here is identical to the code contained in the Data S1 supplement of the paper. Because Github is not designed for file storage and some of the data files exceed Github's upload limits, files in the "Data" folder listed below are not available on Github and can be obtained from the paper supplement here: [Data S1 supplement](https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Feap.1854&file=eap1854-sup-0004-DataS1.zip) (OPENS AS A ZIP FOLDER). Note that you can reproduce the .RDS files in the "Data" folder by running the scripts, but the .RDS files are provided in case you wish to skip this step and dig into the code itself.   

**Content List:**

*Code* [Folder]

  * Appendix-S3-Script.R
  * Simulation-Functions.R
  * Simulation-Script.R
  
*Data* [Folder]

  * *appendix-results* [Folder]
  
      - amdata_object.RDS
      - bias_detection.RDS
      - bias_state.RDS
      - presence_warnings.RDS
      
      
  * *simulation-results* [Folder]
  
      - amdata_object.RDS
      - bias_detection.RDS
      - bias_state.RDS
      - presence_warnings.RDS


**Content Description:** 

*Code*: folder containing three .R files, all three of which are heavily commented to
be followed along with by a user.

  - Appendix-S3-Script.R – an R script to replicate the simulation results
contained in Appendix S3 (copies of which are provided in the Data folder
'appendix-results').
  - Simulation-Functions.R – an R file containing all functions required
to run Appendix-S3-Script.R and Simulation-Script.R.
  - Simulation-Script.R – an R script to replicate the simulation results
contained in the main body of the paper (copies of which are provided in the
Data folder 'simulation-results').

*Data*: folder containing two folders of results

  - *appendix-results*: a folder containing all data produced by the Appendix S3
simulation results (which can be reproduced using Appendix-S3-
Script.R)
    * amdata_object.RDS – an RDS file containing an AMModels
class ‘amData’ object which stores dynamics, encounter histories, data
summaries, and parameter estimates from each of 100 replicates of all
192 appendix simulation scenarios.
    * bias_detection.RDS – an RDS file containing a data.table with
dimensions 19,200 x 18 that stores parameter estimates from each of
100 replicates of all 192 appendix simulation scenarios, for the
occupancy detection parameters p11, p10, and b. This object is used to
produce plots of the bias of detection parameters with the function
simBiasPlot() provided in Simulation-Functions.R and used
in Appendix-Script.R.
    * bias_state.RDS – an RDS file containing a data.table with
dimensions 19,200 x 18 that stores parameter estimates from each of
100 replicates of all 192 appendix simulation scenarios, for the
occupancy state parameters psi (ψ), gamma (γ), and epsilon (ε). This
object is used to produce plots of the bias of state parameters with the
function simBiasPlot() provided in Simulation-Functions.R
and used in Appendix-Script.R.
    * presence_warnings.RDS – an RDS file containing a data.table
storing the scenario name (‘scenario’) and replicate number (‘rep’) of
each scenario-replicate that received a warning during model fitting
from the program PRESENCE. The column ‘conv.warning’ tracks
whether this scenario-replicate received a convergence warning, and if
so, at what value. The ‘VC.warning’ column tracks whether this
scenario-replicate received a warning about the variance-covariance
matrix. This object is used so that scenario-replicates that failed to
converge may be removed from plots of parameter bias with the 
function simBiasPlot() provided in Simulation-Functions.R
and used in Appendix-Script.R.


- *simulation-results*: a folder containing all data produced by the Appendix
S3 simulation results (which can be reproduced using SimulationScript.R)
    * amdata_object.RDS – an RDS file containing an AMModels
class ‘amData’ object which stores dynamics, encounter histories, data
summaries, and parameter estimates from each of 500 replicates of all
128 simulation scenarios.
    * bias_detection.RDS – an RDS file containing a data.table with
dimensions 192,000 x 18 that stores parameter estimates from each of
500 replicates of all 128 simulation scenarios, for the occupancy
detection parameters p11, p10, and b. This object is used to produce
plots of the bias of detection parameters with the function
simBiasPlot() provided in Simulation-Functions.R and used
in Simulation-Script.R.
    * bias_state.RDS – an RDS file containing a data.table with
dimensions 192,000 x 18 that stores parameter estimates from each of
500 replicates of all 128 simulation scenarios, for the occupancy state
parameters psi (ψ), gamma (γ), and epsilon (ε). This object is used to
produce plots of the bias of state parameters with the function
simBiasPlot() provided in Simulation-Functions.R and used
in Simulation-Script.R.
    * presence_warnings.RDS – an RDS file containing a data.table
storing the scenario name (‘scenario’) and replicate number (‘rep’) of
each scenario-replicate that received a warning during model fitting
from the program PRESENCE. The column ‘conv.warning’ tracks
whether this scenario-replicate received a convergence warning, and if
so, at what value. The ‘VC.warning’ column tracks whether this
scenario-replicate received a warning about the variance-covariance
matrix. This object is used so that scenario-replicates that failed to
converge may be removed from plots of parameter bias with the
function simBiasPlot() provided in Simulation-Functions.R
and used in Simulation-Script.R.