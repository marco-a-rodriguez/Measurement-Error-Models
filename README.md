#### Code and data files for "Measurement error models reveal the scale of consumer movements along an isoscape gradient" (Rodríguez 2021)   

#### Data files:   
MEMdataSalmonRasmussen2009.csv   
* _Original source: Rasmussen et al. (2009), Supplementary Material, file JANE_1511_sm_AppendixS1.doc; https://doi.org/10.1111/j.1365-2656.2008.01511.x_  

MEMdataBluesharkBird2018.csv   
* _Original source: Bird et al. (2018); https://datadryad.org/stash/dataset/doi:10.5061/dryad.d1f0d_

MEMdataCatsharkBird2018.csv   
* _Original source: Bird et al. (2018); https://datadryad.org/stash/dataset/doi:10.5061/dryad.d1f0d_ 

#### R code files:
MEM_empirical_data.r
* _R code to generate results of the section "Results: Estimates of movement scale and fractionation in three fish species" in "Measurement error models
reveal the scale of consumer movements along an isoscape gradient"_   
* _Reads data from MEMdataSalmonRasmussen2009.csv, MEMdataBluesharkBird2018.csv, or MEMdataCatsharkBird2018.csv_   
* _Calls Stan programs MEMgaussian.stan, MEMlaplace.stan, or MEMstudent.stan_   
* _Loads required functions from MEMfunctions.r_   

MEM_simulated_data.r   
* _R code to generate results of the section "Results: Performance of estimators in simulated scenarios" in "Measurement error models reveal the scale of consumer movements along an isoscape gradient"_   
* _Calls Stan program MEMsims.stan_   
* _Loads required functions from MEMfunctions.r_   

MEMfunctions.r   
* _R libraries and functions sourced by R programs in "Measurement error models reveal the scale of consumer movements along an isoscape gradient"_   

#### Stan code files:   
MEMgaussian.stan; MEMlaplace.stan; MEMstudent.stan   
* _Stan code to generate results for the section "Applications of the model: quantifying the movements of three fish species" in "Measurement error models reveal the scale of consumer movements along an isoscape gradient"_   
* _Are called by R program MEM_empirical_data.r_   
 
MEMsims.stan   
* _Stan code to generate results of the section "Results: Performance of estimators in simulated scenarios" in "Measurement error models reveal the scale of consumer movements along an isoscape gradient"_   
* _Is called by R program MEM_simulated_data.r_   
