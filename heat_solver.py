import numpy as np
import random
from generate_maze import create_maze
from plot_util import *
import time
import matplotlib.pyplot as plt

# parameters
#seed = int(time.time())
seed = 1668224834 # L = 5
print('seed:', seed)

L = 7
T = 2000
alpha = 1
dx = 1
initial_value = 0.0

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
    max_time = T-1
    t = 0
    inter_coord = [0,0]
    while t < max_time:
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
            
            '''
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
            '''  

        E[t] = energy_sum
        # enforce periodic boundary conditions
        u[t+1,start[0],start[1]] += u[t+1,end[0],end[1]] # generator
        u[t+1,end[0],end[1]] = 0.0 # sink

        t += 1
        if t == T-1:
            print('reached maximum time T =', t)
            break

    return u, E, max_time, inter_coord

def plotoptix_data(u):
    u = u.copy()

    B = np.zeros(u.shape)
    '''
    rise_time = 10
    for t in range(T-1):
        for i,j in idx:
            if u[t+1,i,j] != initial_value and u[t,i,j] == initial_value:
                B[t+1:t+rise_time+1,i,j] = np.linspace(1/rise_time,1,num=rise_time)
                B[t+rise_time+1:t+2*rise_time+1,i,j] = 1 - np.linspace(1/rise_time,1,num=rise_time)
    '''

    tol = 10e-5
    height_change = 0.02
    for t in range(T-1):
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
            diff_count = 0
            if abs(u0 - u1) > tol: diff_count += 1
            if abs(u0 - u2) > tol: diff_count += 1
            if abs(u0 - u3) > tol: diff_count += 1
            if abs(u0 - u4) > tol: diff_count += 1
            
            if  diff_count >= 2:
                if B[t,i,j] < 1:
                    B[t+1,i,j] = B[t,i,j] + height_change
                else:
                    B[t+1,i,j] = B[t,i,j]
            else:
                if B[t,i,j] > 0:
                    B[t+1,i,j] = B[t,i,j] - height_change

    # normalize the values so that the max = 1
    for t in range(T-1):
        for i,j in idx:
            u[t,i,j] = u[t,i,j] * 1 / np.max(u[t])
        u[t,start[0],start[1]] = 1.2
        u[t,end[0],end[1]] = 1.2

    # save B and u in .npy file
    optix_data = np.array([maze, u, B], dtype="object")
    np.save('maze_data_plotoptix.npy', optix_data) 
    
    return B


# solve the PDE
u, E, max_time, inter_coord = SolveHeat(u)


# save data for plotoptix
B = plotoptix_data(u)

# create plot
PlotEnergy(E)
PlotMaze(maze, start, end)
PlotShortestPath2(u[max_time-1], maze, start, end)
#PlotShortestPath1(u[max_time-1], maze, inter_coord, start)
#PlotShortestPath2(u[max_time-1], maze, inter_coord, end)

# create animation
#AnimateSolution(u[:max_time-1], start, end, max_time-1, frame_skip=20, file_name='solution.gif')
AnimateSolution(B, start, end, T-1, frame_skip=20, file_name='change_heat.gif')


