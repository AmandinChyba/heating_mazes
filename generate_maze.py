import numpy
from mazelib import Maze
from mazelib.generate.Prims import Prims
from mazelib.generate.DungeonRooms import DungeonRooms

def create_maze(seed, Lx, Ly):
    seed = seed
    m = Maze(seed)
    #m.generator = Prims(Lx, Ly)
    m.generator = DungeonRooms(Lx, Ly)
    m.generate()
    maze = m.grid
    maze = maze.astype('float32')
    return maze


