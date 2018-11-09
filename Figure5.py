import csv
import matplotlib.pyplot as plt
import math
import numpy as np

plt.figure(5, figsize=(6, 4), dpi=320)

d_array = np.linspace(0,1,100)

betas = [0.8, 0.4, 0.2, 0.1]

plt.plot([], linestyle="None", marker="o", color="#FF00FF", label=r"$\beta = 0.8$")
plt.plot([], linestyle="None", marker="+", color="#000000", label=r"$\beta = 0.4$")
plt.plot([], linestyle="None", marker="x", color="#0000FF", label=r"$\beta = 0.2$")
plt.plot([], linestyle="None", marker="o", markerfacecolor="None", color="#FF0000", label=r"$\beta = 0.1$")

for beta in betas:
    u_norm = [ (beta/(1-beta*beta))*((1/(beta+d*(1-beta)))-(beta+d*(1-beta))) for d in d_array ]

    plt.plot(d_array, u_norm, color="#000000", linestyle="-")

data = []
with open("Figure5_0.8.out", "r") as data_file:
    csvdata_read = csv.reader(data_file, quotechar='"', delimiter=' ')
    for row in csvdata_read:
        if float(row[1]) != 0.:
            data.append([float(row[0]), float(row[1])])

x = [ d[0] for d in data ]
y = [ d[1] for d in data ]

plt.plot(x, y, linestyle="None", marker="o", color="#FF00FF")

data = []
with open("Figure5_0.4.out", "r") as data_file:
    csvdata_read = csv.reader(data_file, quotechar='"', delimiter=' ')
    for row in csvdata_read:
        if float(row[1]) != 0.:
            data.append([float(row[0]), float(row[1])])

x = [ d[0] for d in data ]
y = [ d[1] for d in data ]

plt.plot(x, y, linestyle="None", marker="+", color="#000000")

data = []
with open("Figure5_0.2.out", "r") as data_file:
    csvdata_read = csv.reader(data_file, quotechar='"', delimiter=' ')
    for row in csvdata_read:
        if float(row[1]) != 0.:
            data.append([float(row[0]), float(row[1])])

x = [ d[0] for d in data ]
y = [ d[1] for d in data ]

plt.plot(x, y, linestyle="None", marker="x", color="#0000FF")

data = []
with open("Figure5_0.1.out", "r") as data_file:
    csvdata_read = csv.reader(data_file, quotechar='"', delimiter=' ')
    for row in csvdata_read:
        if float(row[1]) != 0.:
            data.append([float(row[0]), float(row[1])])

x = [ d[0] for d in data ]
y = [ d[1] for d in data ]

plt.plot(x, y, linestyle="None", marker="o", markerfacecolor="None", color="#FF0000")

plt.legend(frameon=False)
plt.xlabel(r"$(r-r_1)/(r_2-r_1)$")
plt.ylabel(r"$u(r)/u_0$")
plt.ylim(0., 1.)
plt.xlim(0., 1.)

plt.savefig("Figure5.png")
