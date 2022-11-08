import numpy as np
import random
from generate_maze import create_maze
from plot_util import *

# parameters
seed = 18
L = 20
T = 1000
alpha = 1
dx = 1

dt = (dx ** 2)/(4 * alpha)
gamma = (alpha * dt) / (dx ** 2)

random.seed(seed)

# create the maze
maze = create_maze(seed,L,L)
idx = np.argwhere(maze != 1)

# randomly generate start/end
start_idx, end_idx = random.sample(range(0, len(idx)), 2)
start = idx[start_idx]
end = idx[end_idx]

# set insulating boundary conditions
u = maze.copy()
u[u == 1] = -1.0
u[u == 0] = 0.5
u = np.tile(u, (T, 1, 1))

# Set the initial conditions
u[0, start[0], start[1]] = 1.0
u[0, end[0], end[1]] = 0.0

def SolveHeat(u):
    E = np.zeros((T-1))
    for k in range(0, T-1):
        energy_sum = 0
        for i,j in idx:
            u0 = u[k][i][j]
            u1 = u[k][i+1][j]
            u2 = u[k][i-1][j]
            u3 = u[k][i][j+1]
            u4 = u[k][i][j-1]
            
            if u1 == -1.0: u1 = u0
            if u2 == -1.0: u2 = u0
            if u3 == -1.0: u3 = u0
            if u4 == -1.0: u4 = u0
            
            u[k+1,i,j] = gamma * (u1 + u2 + u3 + u4 - 4*u0) + u0
            energy_sum += u0

        E[k] = energy_sum
        
        # enforce periodic boundary conditions
        u[k+1,start[0],start[1]] += u[k+1,end[0],end[1]] # generator
        u[k+1,end[0],end[1]] = 0.0 # sink
        #u[k+1,start[0],start[1]] = 1.0 # generator
        #u[k+1,end[0],end[1]] = 0.0 # sink
    
    return u, E

# solve the PDE
u, E = SolveHeat(u)

# create plot
PlotEnergy(E)
PlotMaze(maze, start, end)
PlotShortestPath(u[T-1], maze, start, end)

# create animation
#AnimateSolution(u, start, end, T)


