### Measurement Error Models   
#### Code and data files for "Measurement error models reveal the scale of consumer movements along an isoscape gradient" (Rodríguez 2021)   

##### Data files:   
MEMdataSalmonRasmussen2009.csv   
Original source: Rasmussen et al. (2009), Supplementary Material,   
file JANE_1511_sm_AppendixS1.doc;   
https://doi.org/10.1111/j.1365-2656.2008.01511.x   

MEMdataBluesharkBird2018.csv   
Original source: Bird et al. (2018);   
https://datadryad.org/stash/dataset/doi:10.5061/dryad.d1f0d  

MEMdataCatsharkBird2018.csv   
Original source: Bird et al. (2018);   
https://datadryad.org/stash/dataset/doi:10.5061/dryad.d1f0d  

##### R code files:
MEM_empirical_data.r   
R code to generate results of the section "Results: Estimates of movement scale and fractionation in three fish species" in "Measurement error models
reveal the scale of consumer movements along an isoscape gradient"      
Reads data from MEMdataSalmonRasmussen2009.csv, MEMdataBluesharkBird2018.csv, or MEMdataCatsharkBird2018.csv   
Calls Stan programs MEMgaussian.stan, MEMlaplace.stan, or MEMstudent.stan   
Loads required functions from MEMfunctions.r   

MEM_simulated_data.r   
R code to generate results of the section "Results: Performance of estimators in simulated scenarios" in "Measurement error models reveal the scale of consumer movements along an isoscape gradient"   
Calls Stan program MEMsims.stan   
Loads required functions from MEMfunctions.r   

MEMfunctions.r   
R libraries and functions sourced by R programs in ms "Measurement error models reveal the scale of consumer movements along an isoscape gradient"   

##### Stan code files:   
MEMgaussian.stan   
MEMlaplace.stan   
MEMstudent.stan   
Stan code to generate results for the section "Applications of the model: quantifying the movements of three fish species" in "Measurement error models reveal the scale of consumer movements along an isoscape gradient"   
Are called by R program MEM_empirical_data.r   
 
MEMsims.stan   
Stan code to generate results of the section "Results: Performance of estimators in simulated scenarios" in ms "Measurement error models reveal the scale of consumer movements along an isoscape gradient"   
Is called by R program MEM_simulated_data.r   
