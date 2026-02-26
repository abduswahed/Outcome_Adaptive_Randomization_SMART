# Rcode_GO_SMART
# Updated in April 2024 by Xue Yang
R code to reproduce simulation results in "A Generalized Outcome-Adaptive Sequential Multiple Assignment Randomized Trial Design"

## Contents

- fun.R: This code contains functions for data generation, Monte Carlo replication, and evaluation of the proposed two-stage GO-SMART designs (AR-1 and AR-2), standard SMART design, and RA-SMART design (Wang et al., 2022). Detailed comments are provided for each function.

- simulation_S1.R: This code reproduces the simulation study under scenario S1, including data generation, data analysis, and the figures shown in the manuscript and Web Appendix; 

- three_stage.R: This code is for the three-stage version of the GO-SMART design in comparison with the three-stage SMART design (Web Appendix F). The simulation under the alternative scenario shown in sTable 1 in the Web Appendix can be reproduced by this code.

- SMART_AR_binary.R: This code is for applying the SMART-AR design (Cheung et al., 2015) for binary outcome (Section 5.3.4 and Web Appendix G). The code is modified based on the Rcode provided in the supplementary file of Cheung et al., 2015. The simulation results under S8 scenario when burn-in proportion is 0.5, and Q-function is correctly specified in the Table 2 is reproduced in this code.


## Sample data:

- The "output" folder contains some datasets as examples. They are all generated from the included R code (simulation_S1.R).


## Set up to reproduce the results

The content of this folder should be saved on the user's computer. In all R files, the location be replaced by the location of the folder on the user's computer. The R packages stringr, tidyverse, ggplot2, latex2exp must be installed to reproduce the figures.



R version 4.3.2 (2023-10-31) -- "Eye Holes" was used.






