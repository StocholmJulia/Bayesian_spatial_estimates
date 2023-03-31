clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "", replace text
set more off


/******************************************************************************************
	   Program: lle_est.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Obtain the estimates of the LLE using the 2,5: 50 & 97,5 quantiles bayesian estimates
	   Data:    bc cohort:  cancer.dta"
	            estimates: "winbugs_est.dta"
	   
	   Created: 2022-12-15 by Yuliya Leontyeva 
	   Updated: 2023-03-06 by Yuliya Leontyeva
	   
	   Software: STATA 17
   *****************************************************************************************/  

use "winbugs_est.dta", clear

keep beta* gamma* u* v* 

// create a new frame for the results
frame create temp_est str15(variable) double(q2 q50 q97)

// Save the names of the covariates in the data set
qui ds
local mylist `r(varlist)'


	foreach i of local mylist {
				
		cap egen a25 = pctile(`i'),p(2.5)
		local a25_`i' = a25
		
		cap egen a50 = pctile(`i'),p(50)
		local a50_`i' = a50
		
		
		cap egen a97 = pctile(`i'),p(97.5)
		local a97_`i' = a97
		
		
		frame post temp ("`i'") (`a25_`i'') (`a50_`i'') (`a97_`i'')
		
		drop a25 a50 a97
}

frame change temp_est

frame temp_est: save quantiles_est.dta,replace

// Necessary mata programm:		
mata:

    mata set matastrict off
	
	void gq(string scalar weightsname, string scalar nodesname)
{
	n =  strtoreal(st_local("nodes"))
	i = range(1,n,1)'
	i1 = range(1,n-1,1)'
	muzero = 2
	a = J(1,n,0)
	b = i1:/sqrt(4 :* i1:^2 :- 1)	
	A= diag(a)
	for(j=1;j<=n-1;j++){
		A[j,j+1] = b[j]
		A[j+1,j] = b[j]
	}	
	symeigensystem(A,vec,nodes)
	weights = (vec[1,]:^2:*muzero)'
	weights = weights[order(nodes',1)]
	nodes = nodes'[order(nodes',1)']
	st_matrix(weightsname,weights)
	st_matrix(nodesname,nodes)
}
end


frame rename default estimates
frame estimates {
	use "quantiles_est.dta", clear
	gen ind = _n
	tempfile temp
	save `temp'
}


// Create another frame with cancer cohort:

frame create data 
frame change data
frame data: use cancer.dta, clear

// Merge with popmort file

// Attained age
gen _age = min(floor(age + _t), 89) 

// Attained year 
gen _year = min(floor(year + _t), 2007)

merge m:1 sex _age _year grid using "popmort.dta", nogen keep(match)

// Create splines of ln(_t) and save knots and Rmatrix to use afterwords:

gen lnt = log(_t)
rcsgen lnt, gen(rcs) df(5) orthog 

global ltknots = r(knots)
matrix ltR = r(R)

drop rcs*


// Create splines of age and save knots and Rmatrix to use afterwords:

rcsgen age, df(2) gen(cage) orthog 
global aknots = r(knots)
matrix aR = r(R)

// Copy the existing data set (so that not to destroy) for further manupulations

frame copy data temp 
frame change temp

// keep only necessary variables: 

keep sex year age grid cage* 
sort grid sex year age 
gen ind = _n


// Define all necessary local variables:

local t0 0           // the beginning of the follow-up period
local tmax 70        // 70-year restricted mean survival    
local nodes 30       // # of  intervals follow-up time is divided, usually 30 is enough  
local meanexp exp    // avarege survival for cancer patients if they did not have cancer
local meanobs obs    // average observed survival for cancer patients 
local survprob prob  // name of a covariate in the population life tables 
local maxage 89      // max age in the population life tables
local maxyear 2007   // max year in the population life tables
local grid 478       // # of regions


tempvar S_star_`t0'
gen `S_star_`t0'' =  1  // expected survival at the beginning of follow-up

local using "popmort.dta" // population life tables 

// Loop over each year until tmax 
// Merge on interval specific expected survival at every year
// Calculate cumulative expected survival

    
	local b = `t0' + 1
	
	forvalues i=`b'/`tmax' {		
		
		tempvar S_star_`i'
		
		local j=`i'-1
	
		 gen _age = floor(min(age + `j',`maxage')) 
		 gen _year= floor(min(year  + `j', `maxyear'))  
		
		sort sex grid _year _age

		qui merge m:1 sex grid _year _age using `using', nogen keep(1 3) keepusing(`survprob')

		gen `S_star_`i''=`S_star_`j'' * `survprob'  // cumulative survival 
		
		drop _age _year `survprob' 			
		
		}
	
			
	******************************************************************
	* Find the nodes and the weights for the numerical integration   *
	* using the guassquad command and then calculate the time points *
	******************************************************************

	
    tempname weightsmat nodesmat
	mata gq("`weightsmat'","`nodesmat'") // find optimal abscissas and weights for the adaptive Legendre-Gaussian integration for a specified # of nodes 
	
	// The integral limits must be converted to [-1:1] to apply the Gaussian integration, i.e. we have to calculate:
	// ((b-a)/2)*x_i + (b - a)/2

	forvalues i=1/`nodes' {			// loop over all nodes to create the time points
		local t`i'= (`tmax'-`t0')*0.5*el(`nodesmat',`i',1)+(`tmax'+`t0')*0.5
		tempvar tvar`i'
		 gen `tvar`i''=`t`i''
	}
	
	****************************************************************************
	* Calculate cumulative expected survival at every time point of interest,  *
	* multiply with the weights, and do the integration (summation)        *
	****************************************************************************
		forvalues i=1/`nodes' {
		tempvar S`i' S_W_`i'
		local floort=floor(`t`i'')
		local ceilt=ceil(`t`i'')
		local dist=`t`i''-`floort'
	    gen `S`i''= `S_star_`floort''-(`S_star_`floort''-`S_star_`ceilt'')*`dist' 
		gen `S_W_`i''= `S`i''*el(`weightsmat',`i',1) 
		drop `S`i'' 
	}
	
	/*Sum up to get the mean expected survival*/
	local SW_list `S_W_1'
	forvalues i=2/`nodes' {
		local SW_list `SW_list' `S_W_`i''
	}
	
	 egen `meanexp' = rowtotal(`SW_list') 
	 replace `meanexp' = `meanexp'*(`tmax'-`t0')*0.5 
	
		
// Obtain average of observed survival:
// Loop over each time point and save R(t)S*(t) (the integrand) for each time point       	

// loop over all regions (478)

sort ind 

// merge with the data set with Bayesian estimates 

merge 1:1 ind using `temp', nogen  

forvalues j = 1/`grid' {
		
	tempvar lnt1
	
	gen `lnt1' = log(`tvar1')
	
	qui rcsgen `lnt1' if grid == `j', gen(rcs) knots($ltknots) rmatrix(ltR)
	drop `lnt1'
	
	// loop over all observations in the quantiles data set and create local variables:

        forvalues ind = 1/2 {
				local b`ind'_2 = q2[`ind']
				local b`ind'_50 = q50[`ind']
				local b`ind'_97 = q97[`ind']
			}
				
					   
	    forvalues ind = 3/8 {
			local s = `ind' - 2 
				local gamma`s'_2 = q2[`ind']
				local gamma`s'_50 = q50[`ind']
				local gamma`s'_97 = q97[`ind']		
		} 
		
		
		
		// obtain estimates of random effects:
		        local k = 8 +`j'
				local u`j'_2 = q2[`k']
				local u`j'_50 = q50[`k']
				local u`j'_97 = q97[`k']				
			
		        
				local l = 486 + `j'
				local v`j'_2 = q2[`l']
				local v`j'_50 = q50[`l']
				local v`j'_97 = q97[`l']				
		 	
		
	// obtain log cumulative hazard as a linear predictor at all quantilies at first time point 
	
	foreach q of numlist 2 50 97 {
	tempvar lnH1_`q'
	
	qui gen double `lnH1_`q'' = `gamma1_`q'' + `gamma2_`q''*rcs1 + `gamma3_`q''*rcs2 + `gamma4_`q''*rcs3 + `gamma5_`q''*rcs4 + `gamma6_`q''*rcs5 + `b1_`q''*cage1 + `b2_`q''*cage2 + `u`j'_`q'' + `v`j'_`q''  
	
		
	// obtain relative survival at first time point at different quantilies:
	
	tempvar rs1_`q' predictstat_`q' 
	
	 qui gen double `rs1_`q'' = exp(-exp(`lnH1_`q''))
	 
	 qui gen double `predictstat_`q'' = `rs1_`q''*`S_W_1'
	 
	 drop `lnH1_`q'' `rs1_`q''
	}	
		
		
	drop rcs*	
	forvalues i=2/`nodes' {
		
		tempvar lnt`i'
		 qui gen `lnt`i'' = log(`tvar`i'')
		
		rcsgen `lnt`i'' if grid == `j', gen(rcs) knots($ltknots) rmatrix(ltR) 
		drop `lnt`i''
		
		// also a loop of three quantiles
		
		foreach q of numlist 2 50 97 {
		
		tempvar lnH`i'_`q' 
		
		qui gen double `lnH`i'_`q'' = `gamma1_`q'' + `gamma2_`q''*rcs1 + `gamma3_`q''*rcs2 + `gamma4_`q''*rcs3 + `gamma5_`q''*rcs4 + `gamma6_`q''*rcs5 + `b1_`q''*cage1 + `b2_`q''*cage2 + `u`j'_`q'' + `v`j'_`q''  
		
		tempvar rs`i'_`q'
		qui gen double `rs`i'_`q'' = exp(-exp(`lnH`i'_`q''))
	    
			
		local predictstat_`q' `predictstat_`q'' + (`rs`i'_`q''*`S_W_`i'')
				
		drop `lnH`i'_`q'' 
		 
		}
		
		drop rcs*
	   
	   }
	   
	   	   
	   // Calculate mean observed survival for each quantile
	   
 foreach q of numlist 2 50 97 {	
 qui gen double `meanobs'_`q' = 0.5*`tmax'*(`predictstat_`q'') // approximate integral 
 
 forvalues i=2/`nodes' {
 	drop `rs`i'_`q''
 }
 }

preserve
 
 // save results for each region separately 
keep if grid == `j' 
tempfile temp_`j'
save `temp_`j''

restore 
drop `meanobs'_2 `meanobs'_50 `meanobs'_97 
		
	} 		

// join together all regions 
use `temp_1', clear
	
forvalues j = 2/`grid' {	
append using `temp_`j''  
}


gen lle_2 = exp - obs_2 
gen lle_50 = exp - obs_50
gen lle_97 = exp - obs_97


drop _*

// save data set:

save lle_est.dta, replace
	
log close	


