# Master-Thesis-Public
Public repository for my Master Thesis.

This repository contains all the code as well as output for my master thesis.

The three RDA files contain the results of my analysis. "Final_results_results" contain the main results (so (r)ARMSE and (r)MoaB) This the data which is loaded into the plotexe.rmd file to produce the figures in the thesis. The "Final_results_rawMSE" contains the  1000 MSE values for each estimation method for each scenario. They are not organized but if you want to find the MSEs for a specific scenario then simply check the "Final_results_results" for a scenario number. Then simply multiply that number by 4000, extract the rows from that value - 3999 to that value. "Final_results_Hyperparms" likewise contains the values of the hyper/penalty parameters for each estimator that used them. To find hyperparameters for a specific scenario do the same as with MSEs. 

If you want to replicate my visualizations you only need to run the plot.rmd and it will save every single plot + legend on your computer. 

If you want to re-run the analysis run the analysisexe.rmd It will call all the functions from the "Functions_script.r" and rerun the analysis. The results depend on a set seed so unless you muck that up the results will be the same as what is already avaialble on the github. 

If you want to see how the analysis was conducted simply look into the "Functions_script.r" it contains all the functions from population generation to CV, to estimating ARMSE. Because of the way it was initially written it is probably easiest to read it from bottom up since the end contains the "mother function" which calls all functions together. "analysisexe" is mainly there since it was the easiest way to paralellize the process. it does not do much other than send all the functions into each core and then once completed it organizes them again.

Before the Thesis defense, the goal is that everything should be cleared up and in 3 neat scripts. If you are too early i might not have gotten there yet and it might be messy.

Sorry.
