import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def PlotShortestPath(u_final, maze, start, end):
    u = u_final.copy()
    maze = maze.copy()

    u[u == -1.0] = np.max(u)
    i,j = start
    while i != end[0] or j != end[1]:
        i_old,j_old = i,j

        u1 = [i+1,j]
        u2 = [i-1,j]
        u3 = [i,j+1]
        u4 = [i,j-1]

        if u[u1[0],u1[1]] < u[i,j]: i,j = u1
        if u[u2[0],u2[1]] < u[i,j]: i,j = u2
        if u[u3[0],u3[1]] < u[i,j]: i,j = u3
        if u[u4[0],u4[1]] < u[i,j]: i,j = u4
        
        if i_old == i and j_old == j: 
            break
        
        maze[i,j] = 0.65
    
    maze[maze == 0] = 0.5
    maze[maze == 1] = 0.0
    maze[start[0],start[1]] = 0.3
    maze[end[0],end[1]] = 1.0

    plt.clf()
    plt.pcolormesh(maze, cmap='seismic')
    plt.savefig("maze_solution.png")

def PlotMaze(maze, start, end):
    maze = maze.copy()
    maze[maze == 0] = 0.5
    maze[maze == 1] = 0.0
    maze[start[0],start[1]] = 0.3
    maze[end[0],end[1]] = 1.0
    
    plt.clf()
    plt.pcolormesh(maze, cmap='seismic')
    plt.savefig("maze.png")

def PlotHeatMap(u_k, k, start, end):
    # normalize the values so that the max = 1
    u_k = u_k * 1 / np.max(u_k)
    u_k[start[0],start[1]] = 1.2
    u_k[end[0],end[1]] = 1.2
    
    plt.clf()
    plt.pcolormesh(u_k, cmap='jet', vmin=-0.2)
    plt.colorbar()
    return plt
   
def PlotEnergy(E_k):
    plt.clf()
    plt.title('Total energy over time')
    plt.xlabel('time')
    plt.ylabel('total energy')
    plt.ylim([0.9, 1.1])
    plt.plot(E_k)
    plt.savefig("energy.png")

def AnimateSolution(u, start, end, T, frame_skip=5):
    u = u.copy()

    def animate(k):
        i = k * frame_skip 
        PlotHeatMap(u[i], i, start, end)

    anim = animation.FuncAnimation(plt.figure(), animate, interval=1, 
                                   frames=int(T/frame_skip), repeat=False)
    anim.save("solution.gif")


















