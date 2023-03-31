clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "'", replace text
set more off


/******************************************************************************************
	   Program: adjm.do
	   
	   Project: Examine spatial variation in the loss in life expectancy.  
	   
	   Purpose: Create adjacency matrix 
	   
	   Data:  1259030002_sla06aaust_shape - SLA areas
	    
	   
	   Created: 2023-01-19 by Yuliya Leontyeva  
	   Updated: 2023-02-07 by Yuliya Leontyeva
	            2023-02-22 by Yuliya Leontyeva
				2023-03-08 by Yuliya Leontyeva
	   
	   Software: STATA 17
   *****************************************************************************************/  
// The shape zip file can be downloaded from https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1259.0.30.0022006

// Convert shapefiles to Stata format:

unzipfile 1259030002_sla06aaust_shape.zip, replace 
spshape2dta SLA06aAUST, replace saving(SLA06)
use SLA06, clear
 
// Tell Stata that _CX and _CY are longitude and latitude values 

spset, modify coordsys(latlong, kilometers)

// Create a data set only for Queensland:

keep if STATE_CODE == "3"
// 478 observations

// Explore the data set
*codebook

sort SLA_5DIGIT

// we have to create grid variabel, i.e. a variable with the same name for SLAs as in cancer cohort to be able to match them later:
gen grid = _n

// verify that a new variable grid really does uniquely identify the observations
bysort grid: assert _N ==1

spcompress, force // so that shape file contains the geographical areas only for Queensland

// Tell sp to use the common GRID variable:

spset grid, modify replace 

// Create an adjacency matrix:

spmatrix create contiguity W, replace // our matrix contains 18 islands, i.e. SLAs without neighbours

// Create a variable adj_grid with the number of neigbours for each SLA
spmatrix summarize W, gen(adj_grid)

// Obtain a list of the SLAs without neigbours 
list grid if adj_grid == 0 


timer off 1
timer list 
log close	
