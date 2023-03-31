clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "Z:\Australia\Australia\Study\Log\cif_est_`todaydate'", replace text
set more off



/******************************************************************************************
	   Program: cif_est.do
	   
	   
	   Purpose: Obtain the estimates of the LLE using the 2,5: 50 & 97,5 quantiles bayesian estimates
	   
	   Data:    cancer.dta (cancer cohort)
	            winbugs_est.dta (Estimation result)
			   
			   
	   Created:            by Yuxin Huang
	   Updated: 2023-03-13 by Yuliya Leontyeva	   
	    
	   
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

// mata programm:		
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


	******************************************************************
	*  Import dataset from WinBUGS estimate  *
	*  contains 2%, 50%, 97% of each beta,gamma,u,v*
	******************************************************************

frame rename default estimates
frame estimates {
	use quantiles_est.dta, clear
	gen ind = _n
	
	tempfile temp
	save `temp'
}


// Create another frame (named data) for cancer cohort:


frame create data 
frame change data
frame data: cancer.dta, clear


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

	******************************************************************
	* Step 1: Calculate expected survival S*   *
	* and determine the weight and timepoint for integration *
	******************************************************************

// define all local variables

local t0 0           // the beginning of the follow-up period
local tmax 10        // 10- crude probability of death 
local nodes 30       // # of  intervals follow-up time is divided, usually 30 is enough  
local cr cr          // crude probability of death due to cancer
local netm netm        // net mortality at tmax (1 - relative survival at maximum time)
local survprob prob  // name of a covariate in the population life tables 
local maxage 89      // max age in the population life tables
local maxyear 2007   // max year in the population life tables
local grid 478       // # of regions


tempvar S_star_`t0'
gen `S_star_`t0'' =  1  // expected survival at the beginning of follow-up

local using "`path'\Data\CleanData\popmort.dta"
 // population life tables 

	local b = `t0' + 1
	
	
// Loop 1: this loop is for obtaining cumulative expected survival for each year in 10-year follow-up (S*_0 to S*_10) //
		
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

// 	Find the nodes and the weights for the numerical integration

	tempname weightsmat nodesmat
	mata gq("`weightsmat'","`nodesmat'") // find optimal abscissas and weights for the adaptive Legendre-Gaussian integration for a specified # of nodes 
	
	
// Loop 2: loop over all nodes in Legendre-Gaussian, and generate coorresponding time point
// To be able to change into an integral over [-1,1]   

	forvalues i=1/`nodes' {			// loop over all nodes to create the time points
		local t`i'= (`tmax'-`t0')*0.5*el(`nodesmat',`i',1)+(`tmax'+`t0')*0.5
		tempvar tvar`i'
		 gen `tvar`i''=`t`i''
	}
	
	
// Loop 3: obtain expected survival for each time points(instead of each year), also generate a new varible: S_W_i, which is expected survival times weights at time point (will be used for later integration) 

	*set trace on 
		forvalues i=1/`nodes' {
			tempvar S`i' S_W_`i'
			local floort=floor(`t`i'')
			local ceilt=ceil(`t`i'')
			local dist=`t`i''-`floort'
	    		gen `S`i''= `S_star_`floort''-(`S_star_`floort''-`S_star_`ceilt'')*`dist' 
			gen `S_W_`i''= `S`i''*el(`weightsmat',`i',1) 
			drop `S`i'' 
		
		}
		
	******************************************************************
	* Step 2: Calculate relative survival R   *
	* and obtain lnH, to compute lambda *
	******************************************************************
	
sort ind 

// paste the data set with Bayesian estimates into cancer cohort 

merge 1:1 ind using `temp', nogen

// Save the estimates in local variables 
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


// Obtain estimates for all grids separately

forvalues j = 1/`grid' {
	
	preserve 
	
	// obtain estimates of random effects:
		        local k = 8 +`j'
				local u`j'_2 = q2[`k']
				local u`j'_50 = q50[`k']
				local u`j'_97 = q97[`k']				
			
		        
				local l = 486 + `j'
				local v`j'_2 = q2[`l']
				local v`j'_50 = q50[`l']
				local v`j'_97 = q97[`l']	
	
	
	keep if grid == `j'
		
	tempvar lnt1
	
	gen `lnt1' = log(`tvar1')
	
	qui rcsgen `lnt1', gen(rcs) dgen(drcs) knots($ltknots) rmatrix(ltR)
	drop `lnt1'
	
	// obtain log cumulative hazard as a linear predictor at all quantilies at first time point (log cumulative excess hazard, which is eta in paper)
	
	foreach q of numlist 2 50 97 {
	tempvar lnH1_`q'
	
	
	qui gen double `lnH1_`q'' =  `gamma1_`q'' + `gamma2_`q''*rcs1 + `gamma3_`q''*rcs2 + `gamma4_`q''*rcs3 + `gamma5_`q''*rcs4 + `gamma6_`q''*rcs5 + `b1_`q''*cage1 + `b2_`q''*cage2 + `u`j'_`q'' + `v`j'_`q''  
	

    * derivative of spline function
	tempvar dspline1_`q'
	qui gen `dspline1_`q'' = `gamma2_`q''*drcs1 + `gamma3_`q''*drcs2 + `gamma4_`q''*drcs3 + `gamma5_`q''*drcs4 + `gamma6_`q''*drcs5
	
		
	// obtain relative survival and excess hazard at first time point at different quantilies:
	
	 tempvar rs1_`q' lambda1_`q' crudepc_`q' 
	
	 qui gen double `rs1_`q'' = exp(-exp(`lnH1_`q''))
	
	
  // excess hazard as derivative of d ln(h(lnt)) / dt	
  
	 qui gen double `lambda1_`q'' = 1/`tvar1'*exp(`lnH1_`q'')*`dspline1_`q''
	  
	  
	 qui gen double `crudepc_`q'' = `rs1_`q''*`S_W_1'*`lambda1_`q''
	 
	
	 drop `lnH1_`q'' `lambda1_`q'' `dspline1_`q'' `rs1_`q''
	
	}
	
	drop rcs* drcs*
	
		
	forvalues i=2/`nodes' {
		
		 tempvar lnt`i'
		 qui gen `lnt`i'' = log(`tvar`i'')
		
		 qui rcsgen `lnt`i'', gen(rcs) dgen(drcs) knots($ltknots) rmatrix(ltR) 
		 drop `lnt`i''
		
		// also a loop of three quantiles
		
		foreach q of numlist 2 50 97 {
				
		tempvar lnH`i'_`q' 
		
		qui gen double `lnH`i'_`q'' = `gamma1_`q'' + `gamma2_`q''*rcs1 + `gamma3_`q''*rcs2 + `gamma4_`q''*rcs3 + `gamma5_`q''*rcs4 + `gamma6_`q''*rcs5 + `b1_`q''*cage1 + `b2_`q''*cage2 + `u`j'_`q'' + `v`j'_`q''  
	    
		
		tempvar dspline`i'_`q'
		
		qui gen `dspline`i'_`q'' = `gamma2_`q''*drcs1 + `gamma3_`q''*drcs2 + `gamma4_`q''*drcs3 + `gamma5_`q''*drcs4 + `gamma6_`q''*drcs5
		
		
		tempvar rs`i'_`q' 
		qui gen double `rs`i'_`q'' = exp(-exp(`lnH`i'_`q''))
	    
		// calculate excess hazard at each time point!
		tempvar lambda`i'_`q'
		 
		qui gen double `lambda`i'_`q'' = 1/`tvar`i''*exp(`lnH`i'_`q'')*`dspline`i'_`q''
		
			
		local crudepc_`q' `crudepc_`q'' + (`rs`i'_`q''*`S_W_`i''*`lambda`i'_`q'')
		
				
		drop `lnH`i'_`q'' `dspline`i'_`q'' 
		 
		}
		
		drop rcs* drcs*
	   
	   }
	   
	   	   
	   // Calculate final results of the integration for cr for each quantile
    
 foreach q of numlist 2 50 97 {	
 qui gen double `cr'_`q' = 0.5*`tmax'*(`crudepc_`q'') // approximate integral 
 qui gen double `netm'_`q' =  1 - `rs`nodes'_`q''
 
 forvalues i=2/`nodes' {
 	drop `rs`i'_`q'' `lambda`i'_`q''
 }
 }

 
 // Save the results for each grid in a separate temporary file
tempfile temp_`j'
save `temp_`j''
	  
// Restore the initial data set 	  
restore 
 
 
  }
  
// join together all regions 
use `temp_1', clear
	
forvalues j = 2/`grid' {	
append using `temp_`j''  
}


sort ind 

save cr_est.dta, replace
	
// Turn of the timer and list it:

log close	
