set PGSL; #set of grocery stores (both potential and current, 1,...,64 are current) 
param xp {j in PGSL}; # indices refer to addresses in data document
param yp {j in PGSL}; # latitute and logitude included

set CB; #set of all census block groups
param xc {j in CB}; # indices refer to lat/lon coordinates of centroids
param yc {j in CB};
param pop {j in CB}; # includes population in each census block group

param num_stores >= 0; #determines the number of stores we want to open
param aversion_to_inequality >=0; #aversion to inequality parameter from aversion.mod

var loc_used {PGSL cross CB} binary; 
	# 1 if census block group j is assigned to grocery store i
	
var store_used {j in PGSL} binary;
	#1 if store j is used
	
var w_cb {j in CB}; #to make linear objective function

var x_cb {j in CB}; #x_cb constraint, i.e. sum over stores of distances*location used for each census block group

var d_cb {PGSL cross CB}; #manhattan distance


minimize KPEDE: 1/aversion_to_inequality*log(1/480*(sum {i in CB} w_cb[i]));
	#objective function
	
subject to Assignment {j in CB}: sum {i in PGSL} loc_used[i,j] = 1;
	# ensures each tract is assigned to exactly one store
	
subject to Stores: sum {i in PGSL} store_used[i] <= num_stores;
	#do not exceed the number of stores available
	
subject to w_cb_constraint {j in CB}: w_cb[j]>=  exp(aversion_to_inequality*pop[j]*x_cb[j]);
	#calculates w_cb for linear objective function

subject to x_cb_constraint {j in CB}: x_cb[j] = sum {i in PGSL} loc_used[i,j]*d_cb[i,j];
	#calculates sum over stores of distances*location used for each census block group

subject to d_cb_constraint {i in PGSL, j in CB}: d_cb[i,j] = abs(xp[i] - xc[j]) + abs(yp[i] - yc[j]);
	#calculates Manhattan distance
	
subject to OnlyIfOpen {i in PGSL, j in CB}: loc_used[i,j] <= store_used[i];
	#a store is only used if it is open

subject to First64AreOpen {i in 1..64}: store_used[i] = 1;
	#we use the grocery stores we have currently;

#TO RUN THE PROGRAM:
#option solver BARON;
#model kpede.mod; data data.dat; solve;

#To see what baron is doing:
#option baron_options 'outlev=1'; 

#To make baron run longer:
#option baron_options 'outlev=1 maxtime=2500'; 

