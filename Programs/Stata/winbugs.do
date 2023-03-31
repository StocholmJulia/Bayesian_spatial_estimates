clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "", replace text
set more off


/******************************************************************************************
	   Program: winbugs.do
	   
	   Project: Examine spatial variation in crude probability of death and the loss in life expectancy.  
	   
	   Purpose: Run Bayesin Spatial model within WinBugs
	   
	   Data:    cancer cohort: "cancer.dta"
	            popmort files: "popmort.dta"
	            adjacency matrix: "adj_matrix.dta"
	            model specification: "Model.txt"
	           
	   
	   Created: 2022-12-02 by Yuliya Leontyeva  
	   Updated: 2023-01-13 by Yuliya Leontyeva
	            2023-03-17 by Yuliya Leontyeva 
	   
	   Software: STATA 17
   *****************************************************************************************/  
// Set the timer of the programm

timer on 1   
clear all 

// Preparation of the data for WinBugs 

// rename default frame:

frame rename default work 

//Load prepared cancer data: 

use cancer, clear


// Merge with popmort file

// Attained age
gen _age = min(floor(age + _t), 89) 

// Attained year 
gen _year = min(floor(year + _t), 2007)

// Merge with population file using sex, bagegroup, year & grid and keep only those who are matched :

// This mortality file is obtained after fitting the Poisson model
merge m:1 sex _age _year grid using "popmort.dta", nogen keep(match)

sort grid 
gen ind = _n // This variable will be used to link different data sets

// We have to include extra variables from adjusted matrix data set:

frame create temp 
frame temp {
	use "adj_matrix.dta", clear
	sort grid adj_grid
	gen ind = _n
}

// First, add a variable adj_grid to the existing data set.
// adj_grid variable contains ID numbers of the adjacent areas
// for each area 

// Link two data frames with id:

frlink 1:1 ind, frame(temp)
frget adj_grid, from(temp)

drop temp // drop a variable connecting two frames 

// create a vector of length 478 (the total number of areas) giving # of neighbours for each area:
frame temp {
	sort grid
	by grid: gen num = _N
	by grid: keep if _n == 1
	replace ind = grid 
}

frlink 1:1 ind, frame(temp)

frget num, from(temp)
drop temp 

// Save the data in WinBugs fromat as a list
// WinBugs does not allow variabels with underscore, so we tell WinBugs to rename it

// create splines for age 

rcsgen age, df(2) gen(cage) orthog 

// Fit the model to obtain initial values 
stpm2 cage1 cage2, df(5) scale(hazard) bhazard(rate)

local cage1 `= [xb][cage1]'
local cage2 `= [xb][cage2]'

local rcs1 `= [xb][_rcs1]'
local rcs2 `= [xb][_rcs2]'
local rcs3 `= [xb][_rcs3]'
local rcs4 `= [xb][_rcs4]'
local rcs5 `= [xb][_rcs5]'
local con `= [xb][_cons]'


local se_cage1 `= _se[cage1]'
local se_cage2 `= _se[cage2]'

local se_rcs1 `= _se[_rcs1]'
local se_rcs2 `= _se[_rcs2]'
local se_rcs3 `= _se[_rcs3]'
local se_rcs4 `= _se[_rcs4]'
local se_rcs5 `= _se[_rcs5]'
local se_cons `= _se[_cons]'


cd "WinBugs"

local N = _N


qui wbslist (var grid _d, name(grid d) format(%3.0f)) (var _rcs1 _rcs2 _rcs3 _rcs4 _rcs5 _d_rcs1 _d_rcs2 _d_rcs3 _d_rcs4 _d_rcs5 cage1 cage2 rate, name(rcs1 rcs2 rcs3 rcs4 rcs5 drcs1 drcs2 drcs3 drcs4 drcs5 cage1 cage2 rate) format(%4.3f)) (var _t,name(t) format(%4.3f)) (vector num if num !=.,format(%2.0f)) (vector adj_grid if adj_grid != . , name(adj) format(%2.0f)) (sumNumNeigh = 2724) (Nsla = 478) (datarows = `N') using Data.txt, replace

 
// Create 2 chains of initial values 


frame create temp2
set seed 1368053

frame temp2 {
forvalues i = 1/2 {

// for betas and gamams we take N
    local beta1 = rnormal(`cage1',`se_cage1')
    local beta2 = rnormal(`cage2',`se_cage2')
	
	local gamma1 = rnormal(`con',`se_cons')
	local gamma2 = rnormal(`rcs1',`se_rcs1')
    local gamma3 = rnormal(`rcs2',`se_rcs2')
    local gamma4 = rnormal(`rcs3',`se_rcs3')
    local gamma5 = rnormal(`rcs4',`se_rcs4')
    local gamma6 = rnormal(`rcs5',`se_rcs5')

	set obs 478 
	
	gen u = 0
	gen v = rnormal(0,1)
	
	local tauu 450
    local tauv 450
	
	gen gamma = .
	replace gamma = `gamma1' if _n == 1
	replace gamma = `gamma2' if _n == 2
	replace gamma = `gamma3' if _n == 3
	replace gamma = `gamma4' if _n == 4
	replace gamma = `gamma5' if _n == 5
	replace gamma = `gamma6' if _n == 6
	
	gen beta = .
	replace beta = `beta1' if _n == 1
	replace beta = `beta2' if _n == 2
			
		
	qui wbslist (var u, format(%3.0f)) (var v, format(%4.3f)) (vector gamma if gamma !=.,format(%4.3f)) (vector beta if beta != ., format(%4.3f)) (tauu = `tauu') (tauv = `tauv')  using init`i'.txt,  replace 
	 

	 drop v u gamma beta 
 
}
}


wbsscript using script.txt, replace model("Model.txt") data("Data.txt") initsfile(init1.txt+init2.txt) logfile("log") coda("out") set(gamma beta u v RERsmoothed var_u_marg fracspatial tauu tauv) burnin(50000) updates(100000) thin(5) /*noquit*/ path("`c(pwd)'") dic 

	
wbsrun using script.txt 
wbscoda using out,clear 


// Read data from coda files created by WinBugs into Stata 			 
wbscoda using out, chains(1 2) clear 


// save the data set in .dta file:

save "winbugs_est.dta", replace


// Turn of the timer and list it:

timer off 1
timer list 
log close	




