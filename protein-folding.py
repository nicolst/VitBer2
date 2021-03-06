import numpy as np
import random

import threading

import time
from tabulate import tabulate

import matplotlib.pyplot as plt


plt.ion()


k_b = 1.38064852 * 10**(-23)

def latex_float(f):
    float_str = "{0:.2g}".format(f).lower()
    if "e" in float_str:
        base, exponent = float_str.split("e")
        float_str = r"{0} \cdot 10^{{{1}}}".format(base, int(exponent))
        return float_str
    else:
        return float_str

class Protein():

    def __init__(self, polymer_length, T, d):
        self.polymer_length = polymer_length
        self.T = T
        self.d = d
        self.N = polymer_length + 1 if polymer_length % 2 == 0 else polymer_length
        self.node_pos = self.N // 2
        self.node = (self.node_pos, self.node_pos)

        self.grid = np.zeros((self.N, self.N), dtype=np.int16)
        for i in range(1, polymer_length + 1):
            self.grid[self.node_pos][i-1] = i

        self.zero_mat = np.zeros((self.N, self.N))

        self.U = (np.random.rand(self.N + 1, self.N + 1)*6.93 - 10.4)*10**(-21)
        for i in range(0, self.N + 1):
            for j in range(0, self.N + 1):
                if -1 <= (i - j) <= 1 or i == 0 or j == 0:
                    self.U[i][j] = 0
        #print(self.U)

        self.current_energy = self.calculate_energy()
        self.energies = [self.current_energy]



    def run(self, wait=False, animate=False):
        if animate:
            cm = plt.get_cmap("jet", self.polymer_length)
            cm.set_bad(color='white')
            fig = plt.matshow(np.ma.masked_where(self.grid == 0, self.grid), cmap=cm)
            cb = plt.colorbar(spacing="uniform")
            #labels = np.array([i for i in range(1, self.polymer_length + 1)])
            labels = np.array([1, self.polymer_length])
            loc = (labels + .5) / (self.polymer_length / (self.polymer_length - 1))
            cb.set_ticks(loc)
            cb.set_ticklabels(labels)
            plt.axis('off')
            plt.title(r"$T={3}$K, $d$: $ {0}$ of ${1}$, $E={2}$J".format(0, self.d, latex_float(self.current_energy), self.T))
            if wait:
                plt.show(block=True)
            else:
                plt.draw()
                plt.pause(0.1)
        for i in range(1, self.d + 1):

            #print("D: {0}/{1}".format(i, self.d))

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

                legal_twist = True

                new_grid = rotate_mat + remainder_mat
                new_energy = self.calculate_energy(new_grid)

                if new_energy < self.current_energy:
                    self.grid = new_grid
                    self.current_energy = new_energy

                    #print("-----{1}-----------------{0}-------------------------\n".format(new_energy, np.argwhere(self.grid == 6)))
                    #self.print_structure()
                elif random.random() < np.exp(-(new_energy - self.current_energy)/(k_b * self.T)):
                    self.grid = new_grid
                    self.current_energy = new_energy
                    #print("-------{2}--------------{0}-----------------{1}---------\n".format(new_energy, rnode_value, np.argwhere(self.grid == 6)))
                    #self.print_structure()

            self.energies.append(self.current_energy)

            if animate:
                if wait:
                    fig = plt.matshow(np.ma.masked_where(self.grid == 0, self.grid), cmap=cm)
                    cb = plt.colorbar(spacing="uniform")
                    cb.set_ticks(loc)
                    cb.set_ticklabels(labels)
                    plt.title(r"$T={3}$K, $d$: $ {0}$ of ${1}$, $E={2}$J".format(i, self.d, latex_float(self.current_energy), self.T))
                    plt.axis('off')
                    plt.show(block=True)
                else:
                    fig.set_data(np.ma.masked_where(self.grid == 0, self.grid))
                    plt.title(r"$T={3}$K, $d$: $ {0}$ of ${1}$, $E={2}$J".format(i, self.d, latex_float(self.current_energy), self.T))
                    plt.draw()
                    plt.pause(0.1)
            else:
                print("D: {0}/{1}".format(i, self.d))


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


#protein_1 = Protein(10, 1, 2)
#protein_1.run(wait=True)
#protein_1.print_structure()
#protein_1.print_structure()

#plt.figure(2)
#d = np.arange(0, protein_1.d + 1)
#plt.plot(d, protein_1.energies)

#plt.show(block=True)

protein_1 = Protein(30, 0.000001, 20000)
#protein_1 = Protein(15, 100, 10000)
#protein_1 = Protein(15, 1000, 8000)

protein_1.run()
plt.figure(2)
d = np.arange(0, protein_1.d + 1)
plt.plot(d, protein_1.energies)
plt.title("prot 1")

plt.figure(3)
mean = []
for i in range(0, len(d)):
    slice = protein_1.energies[:i+1]
    mean.append(np.average(slice))
plt.plot(d, mean)

plt.show(block=True)