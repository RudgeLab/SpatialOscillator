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
G1 = nx.DiGraph()
# add cell's ids as nodes
G1.add_node(1)
G1.add_nodes_from(nodes)
# add (parent_id, child_id) tuples as edges
G1.add_edges_from(edges)

# create directed graph
G2 = nx.DiGraph()
# add cell's ids as nodes
G2.add_node(1)
G2.add_nodes_from(nodes)
# add (parent_id, child_id) tuples as edges
G2.add_edges_from(edges)



path = 'SIMS/'
gamma = 0.3
d = 100.0

with open(path+f"data_contSIM_{gamma}_{d}.json") as json_file:  
    database = json.load(json_file)

cell1 = str(4616)
cell2 = str(9855)
fluos1 = []
fluos2 = []
for k in range(392, len(database)):
    print(k)
    
    ## CELL 1
    # hasn't divided
    if cell1 in database[str(k)].keys():
        fluos1.append([k, database[str(k)][cell1]['fluo']])
    # divided
    else:
	    try:
	        succ1 = list(G1.successors(int(cell1)))
	        cell1 = str(succ1[0])
	        fluos1.append([k, database[str(k)][cell1]['fluo']])
	    except:
	        pass
        
    ## CELL 2
    # hasn't divided
    if cell2 in database[str(k)].keys():
        fluos2.append([k, database[str(k)][cell2]['fluo']])
    # divided
    else:
	    try:
                succ2 = list(G2.successors(int(cell2)))
                cell2 = str(succ2[1])
                fluos2.append([k, database[str(k)][cell2]['fluo']])
	    except:
                pass

dic = {'fluo1': fluos1, 'fluo2': fluos2}
with open(f"fluo_track_{gamma}_{d}.json", 'w') as outfile:  
    json.dump(dic, outfile)
