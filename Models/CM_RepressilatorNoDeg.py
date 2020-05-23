import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

#Import Euler integrator for solving ODE system of chemical species inside the cells
from CellModeller.Integration.CLEulerIntegrator import CLEulerIntegrator

max_cells = 150000


def setup(sim):
    # Set biophysics, signalling, and regulation models
    
    #max_sqr
    #biophys = CLBacterium(sim, max_cells=max_cells, jitter_z=False, max_sqs=256**2)
    biophys = CLBacterium(sim, max_cells=max_cells, jitter_z=False)
    integ = CLEulerIntegrator(sim, 3, max_cells)

    # use this file for reg too
    regul = ModuleRegulator(sim)
    # Only biophys and regulation
    sim.init(biophys, regul, None, integ)

    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,0,0)) 

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)

    sim.pickleSteps = 20

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 3.0 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    #cell.growthRate = 0.6
    cell.growthRate = 1.0
    
    # Specify initial concentration of chemical species
    cell.species[:] = [5,0,0]

def specRateCL():
    return '''
    const float d = 198.f;
    const float e = 0.f;
    
    float p1 = species[0];
    float p2 = species[1];
    float p3 = species[2];
    rates[0] = (d + e*p3*p3)/(1 + p3*p3);
    rates[1] = (d + e*p1*p1)/(1 + p1*p1);
    rates[2] = (d + e*p2*p2)/(1 + p2*p2);
    '''

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        #cell.color = [0.0, cell.species[5]/100, 0.0]
        cell.color = [cell.species[0]/500, cell.species[1]/500, cell.species[2]/500]
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.0 + random.uniform(0.0,0.5)
    d2.targetVol = 3.0 + random.uniform(0.0,0.5)

