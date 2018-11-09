import csv
import matplotlib.pyplot as plt
import math
import numpy as np

plt.figure(5, figsize=(6, 4), dpi=320)

d_array = np.linspace(0,1,100)

betas = [0.8, 0.4, 0.2, 0.1]

for beta in betas:
    u_norm = [ (beta/(1-beta*beta))*((1/(beta+d*(1-beta)))-(beta+d*(1-beta))) for d in d_array ]

    plt.plot(d_array, u_norm, color="#000000", linestyle="-")

plt.xlabel(r"$(r-r_1)/(r_2-r_1)$")
plt.ylabel(r"$u(r)/u_0$")
plt.ylim(0., 1.)
plt.xlim(0., 1.)

plt.savefig("Figure5.png")
