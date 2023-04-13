This code produces the data, does the analysis and creates the visualizations that appear on 
Feasibility in MacArthur’s Consumer-Resource Model by  Andrea Aparicio, Tong Wang, Serguei Saavedra, and Yang-Yu Liu.



Folders and files:

	•	Fig1 contains code to create Figure 1

        •       Fig2 contains code to create Figure 2

	•	FeasibilityVolume contains code to produce the data (build synthetic communities and calculate their feasibility volume and other system properties)
	⁃	functions-MACR.R: contains all the necessary functions to run the simulations
	⁃	run-CR-RandomComp_connectance.R: contains routines to build communities governed by the consumer-resource  model and random matrices iteratively, varying their connectance, and using the functions in functions-MACR.R. The output of this script is the feasibility and complexity data used to build figures 3, 4, S3, S7, and S8 
	⁃	run-CR-otherSims.R: contains routines to build communities governed by the consumer-resource  model iteratively, varying their connectance and controlling other system properties (nestedness, different parameter distributions, and magnitude of Z). The The output of this script is the feasibility and complexity data used to build figures 5, S1, S5, and S9

	•	FiguresVolume contains code to do the data analysis and create the visualizations in the manuscript
	⁃	Figure_functions.R: contains all the necessary functions to build the figures
	⁃	Figures_volume.R: reads the data and calls the functions to build the figures





How to generate the data for every result:

1. Larger and more connected communities have a smaller chance of coexistence
	•	Data: run-CR-RandomComp_connectance.R with option: mode = "MCRM_connectance" and mode = "random_competition"

2. The likelihood of feasibility is inversely proportional to the system’s complexity
	•	Data: script: run-CR-RandomComp_connectance.R with option: mode = "MCRM_connectance"
	•	Power law: script Figures_volume.R, function: make_fitted_data

3. The resources growth dynamics determine the complexity-feasibility relationship when n < m
	•	Data: run-CR-RandomComp_connectance.R with option: script: with mode = "n<m"

4. Increased nestedness leads to a loss of feasibility
	•	Data: run-CR-otherSims.R, function: iterate_nest
5. Niche overlap leads to loss of feasibility
	•	Data: run-CR-RandomComp_connectance.R with option: mode = "MCRM_connectance"





