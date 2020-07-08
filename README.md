# LongitSemiCompReproduce
Simulation code and pseudo dataset for reproducibility of results from "Modeling semi-competing risks data as a longitudinal bivariate process" by Nevo et al.

All analyses use the “LongitSemiComp” R package available from https://github.com/daniel258/LongitSemiComp. 

```
install.packages("devtools")
install_github("daniel258/LongitSemiComp")
```

## Data analysis
- In the paper, we analyzed a dataset obtained from the Adult Changes in Thought (ACT) Study, started at 1994. It contains full information on 4,367 people enlisted to the ACT study at the age 65 or later. Upon collecting baseline data, participants have been followed on a biennial basis for updated clinical data and have been getting comprehensive cognitive evaluation for diagnosis of dementia, and specifically, Alzheimer’s disease. 

- We provide a pseudo dataset with similar design with the same sample size and with the same number of covariates, same type of covariates, and same marginal distribution for the covariates as in the ACT data. 

The folder *Data analysis* contains the following files
-	*ACTpseudo.csv*: the pseudo dataset 
-	*ACTpseudo2.5years.R*: All bivariate modeling analyses of semi-competing risks ACT pseudo data under 2.5 year interval partitions. See below for which Tables and Figures the results are analogous to.
-	*ACTpseudo5years.R*: All bivariate modeling analyses of semi-competing risks ACT pseudo data under 5 year interval partitions. See below for which Tables and Figures the results are analogous to.
The results of the analyses are also given in the RData files *ResultsACTpseudo2.5years.RData* and *ResultsACTpseudo5years.RData*

*NOTE* This psuedo data set sole use is to demonstrate how was the ACT data analyzed
The results obtained from analyzing this dataset are not the same as those reported in the paper 

## Simulation Study

There are five type of files associated with the simulation figures (see Instruction to Use):
-	R scripts to run the simulations.
-	RData files storing the raw simulation results
-	R scripts to combine the multiple Rdata files into matrices and data frames
-	Rdata files storing combined matrices and data frames
-	R scripts to reproduce tables in 

The code in the form provided does not need to be run on a cluster; however, for the simulation studies, the code was ran on a cluster. 

-	The subfolder “scripts” contains R scripts to run the simulations, including the seed number. The R script files are of the form “Scen1Cens00N500.R” where “Scen1” indicates the simulation scenario number (1,2 or 3), “Cens” indicates the censoring level (00, 20 or 30) and “N” indicates the sample size (500, 1000, 5000). For N5000, to shorten running time, each 1000 simulation iterations were cut into 10 separate runs of 100 iterations each, with distinct letter after the 5000. E.g. “Scen1Cens00N5000a.R”. These letters are also used in the script to get unique seed number. 
-	Rdata files containing the results of running each script are available from the subfolder “Results”. The naming system of the files is identical to the R script names explained above.
-	Three R scripts (one for each simulation scenario) to combine the multiple files above into 3 csv files, according the 3. These files are named “ResScen1.csv” replacing “Scen1” with “Scen2” or “Scen3” for scenarios 2 or 3.
-	Under the subfolder “Plots and Tables”, each script reproduces a Figure or a Table from the paper. These include Figures A.2 and A.3 and Tables A.4 to A.12.

