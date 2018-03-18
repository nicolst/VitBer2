import numpy as np
import random
import time

import sys, collections
from multiprocessing import Process as Task, Queue

import asyncio

import threading

import time
from tabulate import tabulate

import matplotlib.pyplot as plt


plt.ion()


k_b = 1.38064852 * 10**(-23)



def print_progress(progress):
    sys.stdout.write('\033[2J\033[H') #clear screen
    for filename, percent in progress.items():
        bar = ('=' * int(percent * 20)).ljust(20)
        percent = int(percent * 100)
        sys.stdout.write("%s [%s] %s%%\n" % (filename, bar, percent))
    sys.stdout.flush()


def latex_float(f):
    float_str = "{0:.2g}".format(f).lower()
    if "e" in float_str:
        base, exponent = float_str.split("e")
        float_str = r"{0} \cdot 10^{{{1}}}".format(base, int(exponent))
        return float_str
    else:
        return float_str

class Protein(threading.Thread):

    def __init__(self, polymer_length, T, d, status, name=None):
        super().__init__(name=name)
        self.polymer_length = polymer_length
        self.T = T
        self.d = d
        self.N = polymer_length + 1 if polymer_length % 2 == 0 else polymer_length
        self.node_pos = self.N // 2
        self.node = (self.node_pos, self.node_pos)
        self.status = status

        self.grid = np.zeros((self.N, self.N), dtype=np.int16)
        for i in range(1, polymer_length + 1):
            self.grid[self.node_pos][i-1] = i

        self.zero_mat = np.zeros((self.N, self.N))

        self.U = (np.random.rand(self.N + 1, self.N + 1)*6.93 - 10.4)*10**(-21)
        for i in range(0, self.N + 1):
            for j in range(0, self.N + 1):
                if -1 <= (i - j) <= 1 or i == 0 or j == 0:
                    self.U[i][j] = 0

        self.current_energy = self.calculate_energy()
        self.energies = [self.current_energy]



    def run(self, wait=False):
        for i in range(1, self.d + 1):
            #print("{2} || D: {0}/{1}".format(i, self.d, self.name))
            if self.T == 0.000001:
                self.status.put([self.T, i / self.d])
            legal_twist = False
            while not legal_twist:
                rnode_value = random.randint(2, self.polymer_length - 1)
                rot_dir = 1 if random.randint(0, 1) == 0 else -1

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
                elif random.random() < np.exp(-(new_energy - self.current_energy)/(k_b * self.T)):
                    self.grid = new_grid
                    self.current_energy = new_energy

            self.energies.append(self.current_energy)



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

def d(T, s, d_max):
    return int(d_max * np.exp(-1*s*T))

start = time.time()

s = 1.31 * 10**(-3)
d_max = 30000

temps = np.append(np.array([0.000001]), np.arange(20, 1501, 20))
ds = np.array([d(temp, s, d_max) for temp in temps])

status = Queue()
progress = collections.OrderedDict()

proteins = [Protein(30, temps[i], ds[i], status) for i in range(0, len(temps))]

averages = []
# for i in range(0, 1 + len(proteins) // 4):
#     for j in range(0, 4):
#         if i*4 + j >= len(proteins):
#             break
#         proteins[i*4 + j].start()
#
#     while any(k.is_alive() for k in proteins):
#         time.sleep(1)
#         while not status.empty():
#             protein, percent = status.get()
#             progress[protein] = percent
#             print_progress(progress)
#     for j in range(0, 4):
#         proteins[i*4 + j].join()
#         averages.append(np.average(proteins[i*4+j].energies))

#i = 0
for protein in proteins:
  protein.start()
  #print("Started {0}".format(i))
  #i += 1

while proteins[0].is_alive():
    time.sleep(1)
    while not status.empty():
        protein, percent = status.get()
        progress[protein] = percent
        print_progress(progress)

# protein_1 = Protein(30, 300, 50000, status)
# protein_1.start()
#
#
# while protein_1.is_alive():
#         time.sleep(1)
#         while not status.empty():
#             protein, percent = status.get()
#             progress[protein] = percent
#             print_progress(progress)
#             asyncio.sleep(1)



for protein in proteins:
    protein.join()
    averages.append(np.average(protein.energies))

# protein_1.join()
# space = np.arange(0, 50001)
# for i in range(0, protein_1.d + 1):
#     averages.append(np.average(protein_1.energies[:i+1]))
#
# plt.figure(1)
# plt.plot(space, averages)
#
# plt.figure(2)
# plt.plot(space, protein_1.energies)



end = time.time()

time_taken = end - start

plt.plot(temps, averages)
plt.ylabel(r"$\langle E \rangle$ / J", size=14)
plt.xlabel(r"$T$ / K", size=14)
plt.title("Mean energy versus temperature (took {0}s)".format(time_taken))



print("Took {0} seconds".format(time_taken))

plt.show(block=True)