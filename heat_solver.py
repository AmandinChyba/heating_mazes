import numpy as np
import random
from generate_maze import create_maze
from plot_util import *
import time
import matplotlib.pyplot as plt

# parameters
seed = int(time.time())
#seed = 1667963409
print('seed:', seed)

L = 30
T = 5000
alpha = 1
dx = 1
initial_value = 0.5

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
u[u == 0] = initial_value
u = np.tile(u, (T, 1, 1))

# Set the initial conditions
u[0, start[0], start[1]] = 1.0
u[0, end[0], end[1]] = 0.0

def SolveHeat(u):
    E = np.zeros((T-1))
    
    done = False
    max_time = T
    t = 0
    #u1_spot = []
    #u2_spot = []
    inter_coord = [0,0]
    while t < max_time:
        #u1_spot.append(u[t-2,start[0],start[1]+1])
        #u2_spot.append(u[t-2,start[0],start[1]-1])
        energy_sum = 0
        
        for i,j in idx:
            u0 = u[t][i][j]
            u1 = u[t][i+1][j]
            u2 = u[t][i-1][j]
            u3 = u[t][i][j+1]
            u4 = u[t][i][j-1]
            
            if u1 == -1.0: u1 = u0
            if u2 == -1.0: u2 = u0
            if u3 == -1.0: u3 = u0
            if u4 == -1.0: u4 = u0

            u[t+1,i,j] = gamma * (u1 + u2 + u3 + u4 - 4*u0) + u0
            energy_sum += u0
        
            # check if we should stop running
            if not done:
                smaller = False
                larger = False
                if u1 < initial_value or u2 < initial_value or u3 < initial_value or u4 < initial_value:
                    smaller = True
                if u1 > initial_value or u2 > initial_value or u3 > initial_value or u4 > initial_value:
                    larger = True

                if smaller and larger:
                    max_time = t + 1#(20 - t % 20)
                    done = True
                    inter_coord[0] = i
                    inter_coord[1] = j
                    print('intersection at', i, j)
                    print('max time:', max_time)
                    print("found a path at t =", t)
                

        E[t] = energy_sum
        # enforce periodic boundary conditions
        #u[t+1,start[0],start[1]] += u[t+1,end[0],end[1]] # generator
        u[t+1,start[0],start[1]] = 1.0 # generator
        u[t+1,end[0],end[1]] = 0.0 # sink

        t += 1
        if t == T:
            print('reached maximum time T =', t)
            break

    '''
    plt.plot(u1_spot)
    plt.plot(u2_spot)
    plt.legend(['First line', 'Second line'])
    plt.show()
    
    print('start at', start[0], start[1], u[t-2,start[0],start[1]])
    print(start[0],start[1]+1, u[t-2,start[0],start[1]+1])
    print(start[0],start[1]-1, u[t-2,start[0],start[1]-1])
    '''

    return u, E, max_time, inter_coord

# solve the PDE
u, E, max_time, inter_coord = SolveHeat(u)

# create plot
PlotEnergy(E)
PlotMaze(maze, start, end)
#PlotShortestPath(u[max_time-1], maze, start, end)
PlotShortestPath1(u[max_time-1], maze, inter_coord, start)
PlotShortestPath2(u[max_time-1], maze, inter_coord, end)

# create animation
AnimateSolution(u[:max_time-1], start, end, max_time-1, frame_skip=20)


