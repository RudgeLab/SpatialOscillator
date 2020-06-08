import time
import sys
import os
import math
import json
import pickle
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import CellModeller

# continuous model
def step(p1, p2, p3, gamma, mu, d, e, n, Dt):
    # Update protein concs
    nextp1 = p1 + ((d + e*(p3**n))/(1 + p3**n) - gamma*p1 - mu*p1) * Dt
    nextp2 = p2 + ((d + e*(p1**n))/(1 + p1**n) - gamma*p2 - mu*p2) * Dt
    nextp3 = p3 + ((d + e*(p2**n))/(1 + p2**n) - gamma*p3 - mu*p3) * Dt
    return nextp1, nextp2, nextp3


#### Build lineage's graph from last pickle
# using Deg03 to check correct values of 
path = '../Data/'
model = 'simpleGrowth10_1/'
files = os.listdir(path+model)
files.sort()
# not using module (think on removing it)
files = files[1:]

# get last pickle's lineage
last_pickle = pickle.load(open(path+model+files[-1], 'rb'))

# constructing data
lin = last_pickle['lineage']
nodes = list(lin.keys())
edges = [(v,k) for k,v in lin.items()]

# create directed graph
G = nx.DiGraph()
# add cell's ids as nodes
G.add_node(1)
G.add_nodes_from(nodes)
# add (parent_id, child_id) tuples as edges
G.add_edges_from(edges)


### Loop over pickles

#### Sim params
# Continuous model params
Dt = 0.05
#gamma = 0.3
#d = 198.
e = 0
n = 2
init_conds = [0, 0, 5.0]

# Sim params
pickeSteps = 1
#time between pickles 
p_time = Dt*pickeSteps #hours

gammas = [0.3] #np.logspace(-2, 0, 5, endpoint=True)
ds = [10000.0]#np.logspace(2, 6, 5, endpoint=True)

##############
## MAIN SIM ##
##############
for gamma in gammas:
    for d in ds:
        print(f"Starting sim gamma={gamma}, d={d}")   
        # Initialize pickle 0 with desired initial conditions
        pickle_0 = pickle.load(open(path+model+'step-00000.pickle', 'rb'))
        pickle_0['cellStates'][1].color = init_conds
        pos = pickle_0['cellStates'][1].pos
        pos = [float(i) for i in pos]
        vol = float(pickle_0['cellStates'][1].volume)
        gr = float(pickle_0['cellStates'][1].strainRate / Dt)
	    
	    
        database = {}
        database[0] = {1: {'pos': pos, 'fluo': init_conds}, 'growthRate': gr, 'volume':vol}

        #for idx in range(1, len(files)):
        startproc = time.time()
        for idx in range(1, len(files)):
	    
            #print("######################################################")
            print(f"gamma={gamma}, d={d}, idx: {idx}")
            
	    #if True:
		    #print(f"File SIMS/data_contSIM_{str(gamma)}_{str(d)}.json EXISTS")
		    #continue
		#    print("HOLA")
		#    break
	    #else:
		#    pass
		#print(f"File SIMS/data_contSIM_{str(gamma)}_{str(d)}.json NO EXISTS")
            
            #startwhole = time.time()
            # get pickles lineage
            # this pickle needs to have concs, so get it from pickle save past iteration
            # Create new key in database
            
            database[idx] = {}
            # previous step, stored in database file
            data1 = database[idx-1]
            # current step, stored in pickle file from CM
            data2 = pickle.load(open(path+model+files[idx], 'rb'))

            #print(f"Number of cells: {len(data2['cellStates'].keys())}")
            if len(data2['cellStates'].keys()) > 60000:
                break

            # need to identify which cells divided from one pickle to the next
            # cells present in step 2 but not in 1 (children)
            not_in_prev = list(set(data2['cellStates'].keys()) - set(data1.keys()))

            # cell division
            if len(not_in_prev) > 0:
                # cells in 1 and not in 2 (parents)
                not_in_curr = list(set(data1.keys()) - set(data2['cellStates'].keys()))
                # each cell in not_in_curr (parents)
                for par_cell in not_in_curr:       
                    succ = list(G.successors(par_cell))

                    # create new entries for children in database file
                    database[idx][succ[0]] = {}
                    database[idx][succ[1]] = {}

                    # Update children state based on parent's state
                    # 1) update parent's concentrations from step 1 to 2, to be inherited updated
                    par_conc_old = data1[par_cell]['fluo']
                    p1, p2, p3 = par_conc_old[0], par_conc_old[1], par_conc_old[2]
                    mu = data1[par_cell]['growthRate']
                    nextp1, nextp2, nextp3 = step(p1, p2, p3, gamma, mu, d, e, n, Dt) 

                    # 2) inherit concentration from parent: assign pos, conc, growthRate, volume to each child
                    concs_1 = [nextp1, nextp2, nextp3]
                    new_conc = np.array(concs_1)
                    vol = data1[par_cell]['volume']
                    # Assign pos            
                    # Needs to convert from numpy.float32 to float
                    database[idx][succ[0]]['pos'] = [float(p) for p in data2['cellStates'][succ[0]].pos]
                    database[idx][succ[1]]['pos'] = [float(p) for p in data2['cellStates'][succ[1]].pos]

                    # Assign species conc
                    database[idx][succ[0]]['fluo'] = new_conc.tolist()
                    database[idx][succ[1]]['fluo'] = new_conc.tolist()

                    # Assign growth rate
                    database[idx][succ[0]]['growthRate'] = mu 
                    database[idx][succ[1]]['growthRate'] = mu
                    
                    # Assign volume
                    database[idx][succ[0]]['volume'] = vol/2 
                    database[idx][succ[1]]['volume'] = vol/2

                    # cells remaining in not_in_prev obj
                    new_left = list(set(not_in_prev) - set(succ))

                    # update not_in_prev
                    not_in_prev = new_left

            # no cell division
            else:
                pass

            ##################################
            # RUN CONTINUOUS OR GILLESPIE SIM
            ##################################

            # this is for cells that did't get divided
            #step(p1, p2, p3, gamma, mu, d, e, n, Dt)

            cells_both = list(set(data2['cellStates'].keys()) & set(data1.keys()))
            #start = time.time()
            # FOR EACH CELL
            for c in cells_both:
                database[idx][c] = {}
                # calculate concentrations
                pos_2 = [float(p) for p in data2['cellStates'][c].pos]
                concs = data1[c]['fluo']

                p1, p2, p3 = concs[0], concs[1], concs[2]
                mu = float(data2['cellStates'][c].strainRate / Dt)
                nextp1, nextp2, nextp3 = step(p1, p2, p3, gamma, mu, d, e, n, Dt)
                vol = float(data2['cellStates'][c].volume)
                # assign to next pickle
                database[idx][c]['pos'] = pos_2
                database[idx][c]['fluo'] = [nextp1, nextp2, nextp3]
                database[idx][c]['growthRate'] = mu
                database[idx][c]['volume'] = vol

            #end = time.time()
            #print(f"SIM: {end-start} secs")
            # STORE PICKLE WITH NEW CONCS
                ## Send picke2 (data2) to be read and stored
            # needs to have all concs updated

            #endwhole = time.time()
            #print(f"Time taken: {endwhole-startwhole}")

        endproc = time.time()
        print(f"Whole process gamma={gamma}, d={d} took {endproc-startproc} secs")

        start = time.time()
        print("Storing file")
        with open(f"SIMS/data_contSIM_{str(gamma)}_{str(d)}.json", 'w') as outfile:  
            json.dump(database, outfile)
        end = time.time()
        print(f"saving gamma={gamma}, d={d} took {end-start} secs")
