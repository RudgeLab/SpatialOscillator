import os
import time
import json
import pickle
import numpy as np
import CellModeller

pickle_path = '../Data/simpleGrowth10_1/'
pickle_files = os.listdir(pickle_path)
pickle_files.sort()

sims_path = 'SIMS/'

sims_files = [
    'data_contSIM_0.3_10000.0.json'
]

#idx=str(392) #5000 OK
#idx=str(655) #20000 OK
#idx=str(715) #25000 OK
#idx=str(818) #35000 OK
#idx=str(949) #50000 OK
#idx=str(578)
idx=str(950)
for file in [sims_files[0]]:
    start = time.time()
    with open(sims_path+file, 'rb') as jf:
        database = json.load(jf)
    data2 = pickle.load(open(pickle_path+f"step-00{idx}.pickle", 'rb'))
    end = time.time()
    print(f"Read file took {end-start} secs")
    concR = []
    concG = []
    concB = []

    for k in database[idx].keys():
        concR.append(database[idx][str(k)]['fluo'][0])
        concG.append(database[idx][str(k)]['fluo'][1])    
        concB.append(database[idx][str(k)]['fluo'][2])

    mr = max(concR)
    mg = max(concG)
    mb = max(concB)
    max_RGB = np.array([mr,mg,mb])

    start = time.time()
    for k in database[idx].keys():
        conc = np.array(database[idx][str(k)]['fluo']) / max_RGB
        data2['cellStates'][int(k)].color = conc.tolist()

    with open(file[:-5]+'_'+idx+'.pickle', 'wb') as handle:
        pickle.dump(data2, handle, protocol=pickle.HIGHEST_PROTOCOL)
    end = time.time()
    print(f"Fix and storing file took {end-start} secs")
