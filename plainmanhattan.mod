set PGSL; #set of grocery stores (both potential and current, 1,...,64 are current) 
param xp {j in PGSL}; # indices refer to addresses in data document
param yp {j in PGSL}; # latitute and logitude included

set CB; #set of all census block groups
param xc {j in CB}; # indices refer to lat/lon coordinates of centroids
param yc {j in CB};
param pop {j in CB}; # population of each census block group

param num_stores >= 0; #determines the number of stores we want to open
param aversion_to_inequality >= 0; #aversion to inequality - included here so the same data
								   #file can be used for both programs

var loc_used {PGSL cross CB} binary; 
	#1 if census block group j is assigned to store j
	
var store_used {j in PGSL} binary;
	#1 if store j is opened

var x_cb {j in CB}; #creating x_cb as a shorthand population*distance

var d_cb {PGSL cross CB}; #creating d_cb as a shorthand for the Manhattan distance


minimize Objective: sum {j in CB} pop[j]*x_cb[j]; #minimize total distance traveled to closest grocery store
	
subject to Assignment {j in CB}: sum {i in PGSL} loc_used[i,j] = 1;
	# ensures each census block group is assigned to exactly one store
	
subject to Stores: sum {i in PGSL} store_used[i] <= num_stores;
	# do not exceed the number of stores available

subject to x_cb_constraint {j in CB}: x_cb[j] = sum {i in PGSL} loc_used[i,j]*d_cb[i,j]; 
	# population*distance

subject to d_cb_constraint {i in PGSL, j in CB}: d_cb[i,j] = abs(xp[i] - xc[j]) + abs(yp[i] - yc[j]); 
	#manhattan distance

subject to OnlyIfOpen {i in PGSL, j in CB}: loc_used[i,j] <= store_used[i]; 
	#only use store if it is open

subject to First64AreOpen {i in 1..64}: store_used[i] = 1; 
	#makes sure all existing stores are used in the distribution

	
#TO RUN THE PROGRAM:
#model plainmanhattan.mod; data data.dat; solve;

#TO SEE WHICH STORES ARE USED:
#display store_used;


	