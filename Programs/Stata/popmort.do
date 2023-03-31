clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "", replace text
set more off


/******************************************************************************************
	   Program: popmort.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Create expected mortality rates by borrowing information from the neighbouring areas 
	   
	   Data:    raw mortality file: "pop_dths.txt"
	           
	   
	   Created: 2022-05-01 by Susanna Cramb  
	   Updated: 2022-12-19 by Yuliya Leontyeva
	            2023-01-13 by Yuliya Leontyeva 
				2023-01-23 by Yuliya Leontyeva
				2023-02-09 by Yuliya Leontyeva
	   
	   Software: STATA 17
   *****************************************************************************************/  

// number of areas in whole population mortality file 
local n_grid 478

use "pop_dths", clear 

// In cancer cohort, the population is older 15 years, therefore we can drop agegroups == 1,2,3

drop if agegroup < 4

// Aggregate data by 5-year calendar period 

gen yeargroup = 1 if inrange(year,1997,2002)
replace yeargroup = 2 if inrange(year,2003,2007)

sort yeargroup grid sex agegroup

collapse (sum) pop count, by(grid yeargroup sex agegroup)
tempfile temp
save `temp'


// Read the adjacency matrix 

	use "adj_matrix", clear 
	
	forvalues val = 1/`n_grid' {
	preserve
    qui keep if grid ==`val' // leave only numbers of neighbouring areas for a specific area  
    keep adj_grid
    rename adj_grid grid
    sort grid
	// merge with cancer cohort 
    qui merge 1:m grid using `temp'
	qui keep if _merge==3 | grid==`val' // we keep all observations for a given region and all its neighbours
	 
	 // calculate the total number of deaths and the total population in these areas
	 collapse (sum) pop count, by(sex yeargroup agegroup)
     gen grid=`val'
	 
	 tempfile gridcount`val'
	 qui save `gridcount`val''
	 
	 restore
	}
	 
	 
	 // merge all temporary files:
	
	use `gridcount1', clear
     
forvalues val = 2/`n_grid' {
    append using `gridcount`val''
    }

// sort the data
sort sex grid yeargroup agegroup  

// explore the range of # of deaths and # of population:

codebook count pop 	

// create observed annual mortality rates: 
gen obs_rate= count / pop

// create midage in each age category:
gen midage=agegroup*5 - 3 

codebook midage

// midage from 17-92

replace midage=89 if midage>90

codebook midage

// Create splines of midage with knots at 18 25 50 75 87 89 

rcsgen midage, gen(agercs) knots(18 25 50 75 87 89) 
 
// Fit Poisson model with splines of midage stratified by grid & calendar period

forvalues grid = 1 / `n_grid' {
	forvalues yeargr = 1/2 {
	preserve
	

qui poisson count agercs* if grid == `grid' & yeargroup == `yeargr', exp(pop)
 
// Create the general population mortality rates for each age from 15-89 and calendar year 1997-2007:

clear
range age 15 89 75
if `yeargr' == 1 {
	range year 1997 2002 6
}
else {
	range year 2003 2007 5
}
fillin age year 
drop if age == . | year == .
gen sex = 2
gen grid = `grid'

rcsgen age, gen(agercs) knots(18 25 50 75 87 89)

predict rate, ir

drop agercs*

tempfile rate`grid'`yeargr'
save `rate`grid'`yeargr''

restore
}
}

// Merge all the grids:

use `rate11', clear
append using `rate12'

forvalues grid = 2 / `n_grid' {
	forvalues yeargr = 1/2 {
	append using `rate`grid'`yeargr''
}
}

drop _fillin
sort grid year age 

gen prob = exp(-rate) 

duplicates report grid year age


rename age _age 
rename year _year 
save "popmort", replace 	

log close	