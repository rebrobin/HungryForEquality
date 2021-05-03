set CGSL; # set of current grocery store locations
param xp {j in PGSL}; # indices refer to addresses in data document
param yp {j in PGSL}; # latitute and logitude included

set CB; #set of all census block groups
param xc {j in CB}; # indices refer to lat/lon coordinates of centroids of tracts
param yc {j in CB};
param pop {j in CB}; #include population of each census block group

param num_stores >=0; #number of stores
param beta; #aversion parameter

var aversion; #variable to be calculated

var x_cb {j in CB}; #shorthand for aversion calculation

subject to x_cb_constraint {j in CB}: x_cb[j] = pop[j]*1/num_stores*(sum {i in PGSL} (abs(xp[i] - xc[j]) + abs(yp[i] - yc[j]))); 
	#calculation of x_cb to be used as a shorthand in aversion calculation
		#approximation of the x_cb used in kpede.mod

subject to aversion_constraint: aversion = (sum {j in CB} x_cb[j])/(sum {j in CB} (x_cb[j]*x_cb[j]))*beta; 
	#calculation of aversion to inequality as used in kpede.mod

#TO RUN THE PROGRAM:
#model aversion.mod; data aversion.dat; solve;

#TO SEE THE CALCULATED AVERSION PARAMETER:
#display aversion;