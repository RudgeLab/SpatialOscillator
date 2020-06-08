import json
import os
import time
import numpy as np
import math

path = 'SIMS/'
gammas = [0.5]
ds = [100.0, 1000.0, 100000.0, 1000000.0]
files = []
for gamma in gammas:
    for d in ds:
        files.append(f"data_contSIM_{gamma}_{d}.json")

files_ok = []
files_not_ok = []

idx = 3 #4 total
for file in [files[idx]]:
    print(f"index: {idx}")    
    print(f"File: {file}")
    try:
        files_ok.append(file)
        start = time.time()
        with open(path+file) as json_file:  
            database = json.load(json_file)
        end = time.time()
        print(end-start)
    except:
        print(f"Error with file: {file}")
        files_not_ok.append(file)

        break

    start = time.time()
    database_kymo = []
    width = 5
    for i in np.arange(10, len(database), 10):
        # when used from database file generated in the code, so keys are integers
        #data = database[i]
        # when file is loaded, key are strings
        data = database[str(i)]

        pos_x = np.array([data[k]['pos'][0] for k in data.keys()])
        pos_y = np.array([data[k]['pos'][1] for k in data.keys()])

        """
        xmax = math.ceil(max(list(map(abs, pos_x))))
        ymax = math.ceil(max(list(map(abs, pos_y))))

        grid_size = max(xmax, ymax)
        xx = np.arange(grid_size)
        yy = np.arange(grid_size)
        y,x = np.meshgrid(xx,yy)
        c = grid_size / 2 - 1/2, grid_size / 2-1/2
        r = np.sqrt((x-c[0])**2 + (y-c[1])**2)
        """
        
        
        r_cells = np.sqrt(pos_x**2 + pos_y**2)
        R_cells = np.array([data[k]['fluo'][0] for k in data.keys()])
        G_cells = np.array([data[k]['fluo'][1] for k in data.keys()])
        B_cells = np.array([data[k]['fluo'][2] for k in data.keys()])

#        print(f"Index: {i}, max_rad:{r_cells.max()}")
        nbins = int(r_cells.max() // width)
        bins_acc = []
        R_acc = []
        G_acc = []
        B_acc = []
        for dr in range(nbins):
            bins_acc.append(dr*width)

            idx = np.where((r_cells > dr*width)*(r_cells < (dr+1)*width))
            R_acc.append(np.mean(R_cells[idx]))
            G_acc.append(np.mean(G_cells[idx]))
            B_acc.append(np.mean(B_cells[idx]))

        database_kymo.append([bins_acc, R_acc, G_acc, B_acc])

    end = time.time()
    print(f"creating kymo database took {end-start} secs")
    
    # extracting vectors from file
    x = [d[0] for d in database_kymo]
    R = [d[1] for d in database_kymo]
    G = [d[2] for d in database_kymo]
    B = [d[3] for d in database_kymo]
    
    start = time.time()
    kymo = np.zeros([len(database_kymo), len(x[-1]), 3])
    kymo[:] = np.nan
    for i, x_vals in enumerate(x):
        kymo[i, :len(x_vals), 0] = R[i]
        kymo[i, :len(x_vals), 1] = G[i]
        kymo[i, :len(x_vals), 2] = B[i]

    vals = kymo[:,:,0]
    kymo[:,:,0] = vals / np.nanmax(vals)

    vals = kymo[:,:,1]
    kymo[:,:,1] = vals / np.nanmax(vals)

    vals = kymo[:,:,2]
    kymo[:,:,2] = vals / np.nanmax(vals)

    end = time.time()
    print(f"creating kymo took {end-start} secs")
    
    # save file
    np.save(f"kymos/kymo_{file[13:-5]}.npy", kymo)
