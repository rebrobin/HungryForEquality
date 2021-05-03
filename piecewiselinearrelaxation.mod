set PGSL; #set of grocery stores (both potential and current, 1,...,64 are current) 
param xp {j in PGSL}; # indices refer to addresses in data document
param yp {j in PGSL}; # latitute and logitude included

set CB; #set of all census block groups
param xc {j in CB}; # indices refer to lat/lon coordinates of centroids
param yc {j in CB};
param pop {j in CB}; # includes population in each tract

param num_stores >= 0; #determines the number of stores we want to open
param aversion_to_inequality >=0; #aversion to inequality as computed in aversion.mod

var loc_used {PGSL cross CB} binary; 
	# 1 if CBG j is assigned to grocery store i, 0 otherwise
	
var store_used {j in PGSL} binary;
	#1 if store j is used
	
var w_cb {j in CB}; #for piecewise linear constraints

var x_cb {j in CB}; #x_cb constraint, i.e. sum over stores of distances*location used for each census block group

var b_cb {j in CB}; #maximum of the distances for each census block group centroid

var d_cb {PGSL cross CB}; #manhattan distance
	
minimize KPEDE: sum {i in CB} w_cb[i];
	
subject to Assignment {j in CB}: sum {i in PGSL} loc_used[i,j] = 1;
	# ensures each tract is assigned to exactly one store
	

subject to Stores: sum {i in PGSL} store_used[i] <= num_stores;
	#do not exceed the number of stores available
	
	
#piecewise linear constraints
	
subject to w_cb_0 {j in CB}: w_cb[j]>= aversion_to_inequality*pop[j]*x_cb[j];
	
subject to w_cb_1 {j in CB}: w_cb[j]>=  exp(aversion_to_inequality*pop[j]*b_cb[j]/1000)*(aversion_to_inequality*pop[j]*(x_cb[j]-b_cb[j]/1000)+1); #for x_cb[j] b/w where

subject to w_cb_2 {j in CB}: w_cb[j]>=  exp(aversion_to_inequality*pop[j]*b_cb[j]/100)*(aversion_to_inequality*pop[j]*(x_cb[j]-b_cb[j]/100)+1); #

subject to w_cb_3 {j in CB}: w_cb[j]>=  exp(aversion_to_inequality*pop[j]*b_cb[j]/10)*(aversion_to_inequality*pop[j]*(x_cb[j]-b_cb[j]/10)+1);

#subject to w_cb_4 {j in CB}: w_cb[j]>=  exp(aversion_to_inequality*pop[j]*b_cb[j]/2)*(aversion_to_inequality*pop[j]*(x_cb[j]-b_cb[j]/2)+1);

subject to w_cb_br {j in CB}: w_cb[j]>=  exp(aversion_to_inequality*pop[j]*b_cb[j])*(aversion_to_inequality*pop[j]*(x_cb[j]-b_cb[j])+1);


#normal other constraints (as in kpede.mod)

subject to x_cb_constraint {j in CB}: x_cb[j] = sum {i in PGSL} loc_used[i,j]*d_cb[i,j];

subject to d_cb_constraint {i in PGSL, j in CB}: d_cb[i,j] = abs(xp[i] - xc[j]) + abs(yp[i] - yc[j]);

subject to OnlyIfOpen {i in PGSL, j in CB}: loc_used[i,j] <= store_used[i];

subject to First64AreOpen {i in 1..64}: store_used[i] = 1;

#needed to compute the max of the distances for each resident in a residential area

subject to b_cb_constraint {j in CB}: b_cb[j] = max{i in PGSL} d_cb[i,j]; 


#TO RUN THE PROGRAM:
#option solver BARON;
#model piecewiselinearrelaxation.mod; data data.dat; solve;

#To see what baron is doing:
#option baron_options 'outlev=1'; 

#To make baron run longer:
#option baron_options 'outlev=1 maxtime=2500'; 


