clear all 
/******************************************************************************************
	   Program: master.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Create a master file with all steps 
	           
	   Created: 2022-12-20 by Yuliya Leontyeva  
	   Updated: 2023-01-13 by Yuliya Leontyeva
	            2023-02-27 by Yuliya Leontyeva 
	   
	   Software: STATA 17
   *****************************************************************************************/  

// Set the working directory

cd ""

// Step 1: Load a shape file for Australia, create shape file only for Queensland, create weigning matrix, which will be used to create SLA and a list of neigbours for each SLA 

do adjm.do


// Step 2: Create general population mortality files. The number of death and the total population is calculated by borrowing information from the neighbouring areas

do popmort.do


// Step 4 Obtaining estimates with WunBugs 

do winbugs.do

   
// Step 6 Obtaining estimates of LLE & Cr

do lle_est.do

do cr_est.do
