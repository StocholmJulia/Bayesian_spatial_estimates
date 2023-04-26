To run WinBugs within Stata, please follow the instructions here:
https://www.stata-press.com/data/bas.html


Name of the project: 

Geographical patterns of the loss in life expectancy and 
crude probability of death for cancer patients
 
Description: This is a collaborative work of Yuliya Leontyeva and Yuxin Huang. 
             The purpose of the project is to estimate the loss in life expectancy and
             crude probability of death using full Bayesian flexible parametric relative
             survival spatial model.
             The paper uses a publicly available simulated datasets
             
Date: 2023-03-31 

Created by: Yuliya Leontyeva, PhD student at MEB

Questions / comments regarding Stata code should be addressed yo Yuliya Leontyeva
                     email: yuliya.leontyeva@ki.se
Questions / comments regarding R code should be addressed Yuxin Huang
                     email: yuxin.huang@hdr.qut.edu.au

--------------------------------------------------------------
This document provides explanations regarding the folders included 
on GitHub

--------------------------------------------------------------
FOLDERS:
1. Data

This folder includes:

- cancer.dta  -  contains simulated data to represent 
              a hypothetical gender-specific cancer

Data columns: cancer.txt (15051 records)
====================================================================================
id		Unique identifier for each case.
sex		All values are 2 (one gender).
year		Year of diagnosis. Ranges from 1997 to 2007.
site10group	All values are 1 (hypothetical cancer).
dx_date		Date of diagnosis (DDMMMYYYY).
dth_date	Date of death (DDMMMYYYY).
age		Age at diagnosis (years).
bagegroup	Broad age group at diagnosis (1=0-49 years, 2=50-69 years, 3=70-89 years).
		Missing means age was 90+ years.
grid		Represents the geographical area of residence at diagnosis (values 1 to 478). 
fu_date		Date of censoring (31Dec2007) or date of death (DDMMMYYYY).
exit		Same as fu_date (DDMMMYYYY).


- pop_dths.txt - contains simulated data to represent hypothetical 
                 population mortality files

Data columns: pop_dths.txt (99902 records)
=====================================================================================
year		Ranges from 1997 to 2007.
sex		All values are 2 (one gender).
agegroup	Values from 1 to 19 representing 5-year age groups (0-4,5-9...,90+).
grid		Represents the geographical area of residence at diagnosis (values 1 to 478). 
count		Number of deaths (all causes).
pop		Population size.



- model.txt  -  a file, which contains the likelihood of the FPRM used
                in the analysis. It contains also all prior distributions
                for all parameters in the model



- 1259030002_sla06aaust_shape - zip file dowloaded from Australian Statistiska bureur (archive)

 Data columns: (478 records)
====================================================================================
_ID_orig	    Unique identifier for SLA
_CX_orig            x-coordinate of SLA
_CY_orig            y-coordinate of SLA
STATE_CODE_orig     State code. For all recordes it is 3
SLA_CODE06_orig     SLA_CODE
SLA_NAME06_orig     SLA_NAME
SLA_5DIGIT_orig     5-digit code for SLA
grid_orig           Order number of SLA from 1-478
adj_grid_orig       The number of neigbours for each SLA 


- adj_matrix -      contains the index of SLA and adjacent areas for each SLA. 
                    18 SLAs without neigbours were assigned to the nearest SLAs by
                    geographical proximity and using relevant locality-specific 
                    information. Number of observations =  2724 




2. Programs 

This folder includes Stata & R programs

In Stata folder there are the following do files:

- set_winbugs - Set WinBugs on your computer 

- master.do - the master file where all the below-mentioned do files are run in the order

- adjm - Create weignening matrix W based on the common border and corner

- popmort.do - creates population mortality rates file in Stata format, where the information from 
               SLAs neigbours is used 

- winbugs - runs Bayesian analysis in WinBugs and saves the estimates

- lle_est - obtains estimates for the Loss in Life Expectancy 

- cr_est - obtains estimates for the crude mortality   


 


  
