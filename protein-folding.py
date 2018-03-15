import numpy as np
import random

import time
from tabulate import tabulate
import matplotlib.pyplot as plt

plt.ion()


k_b = 1.38064852 * 10**(-23)

class Protein:

    def __init__(self, polymer_length, T, d):
        self.polymer_length = polymer_length
        self.T = T
        self.d = d
        self.N = polymer_length + 3 if polymer_length % 2 == 0 else polymer_length + 2
        self.node_pos = self.N // 2
        self.node = (self.node_pos, self.node_pos)

        self.grid = np.zeros((self.N, self.N), dtype=np.int16)
        for i in range(1, polymer_length + 1):
            self.grid[self.node_pos][i] = i

        self.zero_mat = np.zeros((self.N, self.N))

        self.U = (np.random.rand(self.N + 1, self.N + 1)*6.93 - 10.4)*10**(-21)
        for i in range(0, self.N + 1):
            for j in range(0, self.N + 1):
                if -1 <= (i - j) <= 1 or i == 0 or j == 0:
                    self.U[i][j] = 0
        print(self.U)

        self.current_energy = self.calculate_energy()



    def run(self):
        cm = plt.cm.viridis
        cm.set_bad(color='white')
        fig = plt.matshow(np.ma.masked_where(self.grid == 0, self.grid), cmap=cm)
        plt.draw()
        plt.pause(0.1)
        for i in range(0, self.d):
            fig.set_data(np.ma.masked_where(self.grid == 0, self.grid))
            plt.title("d: {0}/{1}".format(i, self.d))
            plt.draw()
            plt.pause(0.1)
            print("D: {0}/{1}".format(i, self.d))

            legal_twist = False
            while not legal_twist:
                #if move_count == (self.polymer_length - 2) * 2:
                #    move_count = -1
                #    print("bah")
                #    break
                rnode_value = random.randint(2, self.polymer_length - 1)
                rot_dir = 1 if random.randint(0, 1) == 0 else -1

                #if rnode_value in moves_done.keys():
                #    if len(moves_done[rnode_value]) == 2:
                #        continue
                #    elif rot_dir in moves_done[rnode_value]:
                #       rot_dir *= -1
                #        moves_done[rnode_value].append(rot_dir)
                #else:
                #    moves_done[rnode_value] = [rot_dir]

                rnode_coord_arr = np.argwhere(self.grid == rnode_value)[0]
                rnode_coord = [rnode_coord_arr[0], rnode_coord_arr[1]]

                rotate_mat = self.grid.copy()

                if rnode_value <= self.grid[self.node_pos, self.node_pos]:
                    rotate_mat[rotate_mat > rnode_value] = 0
                else:
                    rotate_mat[rotate_mat < rnode_value] = 0

                remainder_mat = self.grid - rotate_mat

                trans_dist = (self.node_pos - rnode_coord[0], self.node_pos - rnode_coord[1])

                rotate_mat = self.translate(rotate_mat, trans_dist[0], trans_dist[1])
                rotate_mat = np.rot90(rotate_mat, rot_dir, (0, 1))
                rotate_mat = self.translate(rotate_mat, -trans_dist[0], -trans_dist[1])

                combined = np.multiply(remainder_mat, rotate_mat)

                # If there are any non-zero elements in combined, this is not a valid twist
                if combined.any():
                    continue

                new_grid = rotate_mat + remainder_mat
                new_energy = self.calculate_energy(new_grid)

                if new_energy < self.current_energy:
                    legal_twist = True
                    self.grid = new_grid
                    self.current_energy = new_energy

                    print("----------------------{0}-------------------------\n".format(new_energy))
                    self.print_structure()
                elif random.random() < np.exp(-(new_energy - self.current_energy)/(k_b * self.T)):
                    legal_twist = True
                    self.grid = new_grid
                    self.current_energy = new_energy

                    print("---------------------{0}-----------------{1}---------\n".format(new_energy, rnode_value))
                    self.print_structure()
            #if move_count == -1:
            #    return i

    def print_structure(self, grid=None):
        if grid is None:
            grid = self.grid
        printable = grid.tolist()
        for i in range(0, len(printable)):
            for j in range(0, len(printable[0])):
                if printable[i][j] == 0:
                    printable[i][j] = None  #
        print(tabulate(printable, tablefmt="plain", missingval=" "))

    def translate(self, m, x_shift, y_shift):
        return np.roll(m, (x_shift, y_shift), (0, 1))

    def test_rot(self):
        self.grid = self.translate(self.grid, 2, 2)

    def nearest_neighbours(self, x, y, grid=None):
        if grid is None:
            grid = self.grid
        E = 0
        if x + 1 < len(grid[0]):
            E += self.U[(grid[y][x + 1])][(grid[y][x])]
        if x - 1 >= 0:
            E += self.U[grid[y][x - 1]][grid[y][x]]
        if y + 1 < len(grid):
            E += self.U[grid[y + 1][x]][grid[y][x]]
        if y - 1 >= 0:
            E += self.U[grid[y - 1][x]][grid[y][x]]
        return E

    def calculate_energy(self, grid=None):
        if grid is None:
            grid = self.grid
        E = 0
        for i in range(1, self.polymer_length + 1):
            pos = np.argwhere(self.grid == i)[0]
            E += 0.5 * self.nearest_neighbours(pos[1], pos[0], grid)
        return E


protein_1 = Protein(150, 373, 500)
print(protein_1.run())
#protein_1.print_structure()
protein_1.print_structure()