import time
import sys
import os
import math
import random
import json
import pickle
import numpy as np
import pandas as pd
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


path = 'SIMS/'
gamma = 0.3
d = 10000.0
with open(path+f"data_contSIM_{gamma}_{d}BK.json") as json_file:  
    database = json.load(json_file)

data = database['392']
idx = [k for k in data.keys()]
cells = np.random.choice(idx, 1000, replace=False)

cols = list(range(0, len(cells)))
cols.append('rad_edge')

start = 392
df = pd.DataFrame(index=range(start, 1026), columns=cols)
for it in range(start, 1026):
    print("####################################")
    print(f"it: {it}")
    data = database[str(it)]
    posx = np.array([data[ky]['pos'][0] for ky in data.keys()])
    posy = np.array([data[ky]['pos'][1] for ky in data.keys()])
    r_cells = np.sqrt(posx**2+posy**2)
    r_max = r_cells.max()
    
    print(f"Rad max: {r_max}")
    
    df.loc[it]['rad_edge'] = r_max
    # for all cells in "cells"
    for i, cell in enumerate(cells):
        # it hasn't divided
        #print(f"Cell {cell}")
        if cell in data.keys():            
            fluo = data[cell]['fluo']
            pos = data[cell]['pos']
            gr = data[cell]['growthRate']
            #vol = data[cell]['volume']
            rad = np.sqrt(pos[0]**2+pos[1]**2)
            df.loc[it][i] = {'rad':rad, 'fluo': fluo, 'gr': gr}#, 'vol':vol}
            #print(f"Rad cell {cell}: {rad}")
        # divided
        else:
            succ = list(G1.successors(int(cell)))
            #index = np.where(cells==cell)
            #print(f"Cell antes {cells[index]}")
            new_cell = str(succ[random.choice([0,1])])
            cells[i] = new_cell
            #print(f"Cell después {cells[index]}")
            cell = new_cell

            fluo = data[cell]['fluo']
            pos = data[cell]['pos']
            gr = data[cell]['growthRate']
            #vol = data[cell]['volume']
            rad = np.sqrt(pos[0]**2+pos[1]**2)
            df.loc[it][i] = {'rad':rad, 'fluo': fluo, 'gr': gr}#, 'vol':vol}
            #print(f"Rad cell {cell}: {rad}")

df.to_json('fluo_tracking_1000.json')
