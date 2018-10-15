import csv
import matplotlib.pyplot as plt
import math


plt.figure(3, figsize=(12,12), dpi=160)

# Load Figure 3 nonconvex data
nonconvex_data = []
first = True
with open("Figure3_nonconvex.dat", "r") as nonconvex_datafile:
    csvdata_read = csv.reader(nonconvex_datafile, quotechar='"', delimiter=',')
    for row in csvdata_read:
        if first:
            first = False
        else:
            nonconvex_data.append([float(row[0]), float(row[2]), float(row[3]), float(row[4])])

# Arrange data by gamma/tau.

nonconvex_data_dict = {}

# Add data
for row in nonconvex_data:
    gamma = row[1]
    if gamma not in nonconvex_data_dict:
        nonconvex_data_dict[gamma] = {}
    tau = row[2]
    if tau not in nonconvex_data_dict[gamma]:
        nonconvex_data_dict[gamma][tau] = {'h':[], 'e':[]}
    h = row[0]
    e = row[3]
    if h not in nonconvex_data_dict[gamma][tau]['h']:
        nonconvex_data_dict[gamma][tau]['h'].append(h)
        nonconvex_data_dict[gamma][tau]['e'].append(e)

# Sort data
for gamma in nonconvex_data_dict:
    for tau in nonconvex_data_dict[gamma]:
        hlist = nonconvex_data_dict[gamma][tau]['h']
        elist = nonconvex_data_dict[gamma][tau]['e']
        # From https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
        hsorted, esorted = (list(t) for t in zip(*sorted(zip(hlist, elist))))
        nonconvex_data_dict[gamma][tau]['h'] = hsorted
        nonconvex_data_dict[gamma][tau]['e'] = esorted

# Load Figure 3 convex data
convex_data = []
first = True
with open("Figure3_convex.dat", "r") as convex_datafile:
    csvdata_read = csv.reader(convex_datafile, quotechar='"', delimiter=',')
    for row in csvdata_read:
        if first:
            first = False
        else:
            convex_data.append([float(row[0]), float(row[2]), float(row[3]), float(row[4])])

# Arrange data by gamma/tau.

convex_data_dict = {}

# Add data
for row in convex_data:
    gamma = row[1]
    if gamma not in convex_data_dict:
        convex_data_dict[gamma] = {}
    tau = row[2]
    if tau not in convex_data_dict[gamma]:
        convex_data_dict[gamma][tau] = {'h':[], 'e':[]}
    h = row[0]
    e = row[3]
    if h not in convex_data_dict[gamma][tau]['h']:
        convex_data_dict[gamma][tau]['h'].append(h)
        convex_data_dict[gamma][tau]['e'].append(e)

# Sort data
for gamma in convex_data_dict:
    for tau in convex_data_dict[gamma]:
        hlist = convex_data_dict[gamma][tau]['h']
        elist = convex_data_dict[gamma][tau]['e']
        # From https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
        hsorted, esorted = (list(t) for t in zip(*sorted(zip(hlist, elist))))
        convex_data_dict[gamma][tau]['h'] = hsorted
        convex_data_dict[gamma][tau]['e'] = esorted

plt.subplot(321)

plt.plot([], color="#FF0000", linestyle="-", marker="x", label="$\\tau=0.505$")
plt.plot([], color="#02FF02", linestyle="-", marker="o", markerfacecolor="None", label="$\\tau=0.6$")
plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")
plt.plot([], color="#FF75FF", linestyle="-", marker="s", markerfacecolor="None", label="$\\tau=3$")

plt.loglog(nonconvex_data_dict[0.25][0.505]['h'], nonconvex_data_dict[0.25][0.505]['e'], basex=10, basey=10, color="#FF0000", linestyle="-", marker="x")
plt.loglog(nonconvex_data_dict[0.25][0.6]['h'], nonconvex_data_dict[0.25][0.6]['e'], basex=10, basey=10, color="#02FF02", linestyle="-", marker="o", markerfacecolor="None")
plt.loglog(nonconvex_data_dict[0.25][1.0]['h'], nonconvex_data_dict[0.25][1.0]['e'], basex=10, basey=10, color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.loglog(nonconvex_data_dict[0.25][2.0]['h'], nonconvex_data_dict[0.25][2.0]['e'], basex=10, basey=10, color="#000000", linestyle="-", marker="v", markerfacecolor="None")
plt.loglog(nonconvex_data_dict[0.25][3.0]['h'], nonconvex_data_dict[0.25][3.0]['e'], basex=10, basey=10, color="#FF75FF", linestyle="-", marker="v", markerfacecolor="None")

plt.legend()
plt.xlabel(r"$log_{10}h$")
plt.ylabel(r"$log_{10}E_r$")
plt.ylim(1e-4,4e-1)
plt.xlim(1e-2, 3e-1)

plt.subplot(322)

plt.plot([], color="#FF0000", linestyle="-", marker="x", label="$\\tau=0.505$")
plt.plot([], color="#02FF02", linestyle="-", marker="o", markerfacecolor="None", label="$\\tau=0.6$")
plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")
plt.plot([], color="#FF75FF", linestyle="-", marker="s", markerfacecolor="None", label="$\\tau=3$")

plt.loglog(convex_data_dict[0.25][0.505]['h'], convex_data_dict[0.25][0.505]['e'], basex=10, basey=10, color="#FF0000", linestyle="-", marker="x")
plt.loglog(convex_data_dict[0.25][0.6]['h'], convex_data_dict[0.25][0.6]['e'], basex=10, basey=10, color="#02FF02", linestyle="-", marker="o", markerfacecolor="None")
plt.loglog(convex_data_dict[0.25][1.0]['h'], convex_data_dict[0.25][1.0]['e'], basex=10, basey=10, color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.loglog(convex_data_dict[0.25][2.0]['h'], convex_data_dict[0.25][2.0]['e'], basex=10, basey=10, color="#000000", linestyle="-", marker="v", markerfacecolor="None")
plt.loglog(convex_data_dict[0.25][3.0]['h'], convex_data_dict[0.25][3.0]['e'], basex=10, basey=10, color="#FF75FF", linestyle="-", marker="v", markerfacecolor="None")

plt.legend()
plt.xlabel(r"$log_{10}h$")
plt.ylabel(r"$log_{10}E_r$")
plt.ylim(1e-4,4e-1)
plt.xlim(1e-2, 3e-1)


plt.subplot(323)

plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")

plt.loglog(nonconvex_data_dict[0.75][1.0]['h'], nonconvex_data_dict[0.75][1.0]['e'], basex=10, basey=10, color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.loglog(nonconvex_data_dict[0.75][2.0]['h'], nonconvex_data_dict[0.75][2.0]['e'], basex=10, basey=10, color="#000000", linestyle="-", marker="v", markerfacecolor="None")

plt.legend()
plt.xlabel(r"$log_{10}h$")
plt.ylabel(r"$log_{10}E_r$")
plt.ylim(1e-4,4e-1)
plt.xlim(1e-2, 3e-1)

plt.subplot(324)

plt.plot([], color="#FF0000", linestyle="-", marker="x", label="$\\tau=0.505$")
plt.plot([], color="#02FF02", linestyle="-", marker="o", markerfacecolor="None", label="$\\tau=0.6$")
plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")
plt.plot([], color="#FF75FF", linestyle="-", marker="s", markerfacecolor="None", label="$\\tau=3$")

plt.loglog(convex_data_dict[0.75][0.505]['h'], convex_data_dict[0.75][0.505]['e'], basex=10, basey=10, color="#FF0000", linestyle="-", marker="x")
plt.loglog(convex_data_dict[0.75][0.6]['h'], convex_data_dict[0.75][0.6]['e'], basex=10, basey=10, color="#02FF02", linestyle="-", marker="o", markerfacecolor="None")
plt.loglog(convex_data_dict[0.75][1.0]['h'], convex_data_dict[0.75][1.0]['e'], basex=10, basey=10, color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.loglog(convex_data_dict[0.75][2.0]['h'], convex_data_dict[0.75][2.0]['e'], basex=10, basey=10, color="#000000", linestyle="-", marker="v", markerfacecolor="None")
plt.loglog(convex_data_dict[0.75][3.0]['h'], convex_data_dict[0.75][3.0]['e'], basex=10, basey=10, color="#FF75FF", linestyle="-", marker="v", markerfacecolor="None")

plt.legend()
plt.xlabel(r"$log_{10}h$")
plt.ylabel(r"$log_{10}E_r$")
plt.ylim(1e-4,4e-1)
plt.xlim(1e-2, 3e-1)


plt.subplot(325)

plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")

plt.loglog(nonconvex_data_dict[1.0][1.0]['h'], nonconvex_data_dict[1.0][1.0]['e'], basex=10, basey=10, color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")

plt.legend()
plt.xlabel(r"$log_{10}h$")
plt.ylabel(r"$log_{10}E_r$")
plt.ylim(1e-4,4e-1)
plt.xlim(1e-2, 3e-1)

plt.subplot(326)

plt.plot([], color="#FF0000", linestyle="-", marker="x", label="$\\tau=0.505$")
plt.plot([], color="#02FF02", linestyle="-", marker="o", markerfacecolor="None", label="$\\tau=0.6$")
plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")
plt.plot([], color="#FF75FF", linestyle="-", marker="s", markerfacecolor="None", label="$\\tau=3$")

plt.loglog(convex_data_dict[1.0][0.505]['h'], convex_data_dict[1.0][0.505]['e'], basex=10, basey=10, color="#FF0000", linestyle="-", marker="x")
plt.loglog(convex_data_dict[1.0][0.6]['h'], convex_data_dict[1.0][0.6]['e'], basex=10, basey=10, color="#02FF02", linestyle="-", marker="o", markerfacecolor="None")
plt.loglog(convex_data_dict[1.0][1.0]['h'], convex_data_dict[1.0][1.0]['e'], basex=10, basey=10, color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.loglog(convex_data_dict[1.0][2.0]['h'], convex_data_dict[1.0][2.0]['e'], basex=10, basey=10, color="#000000", linestyle="-", marker="v", markerfacecolor="None")
plt.loglog(convex_data_dict[1.0][3.0]['h'], convex_data_dict[1.0][3.0]['e'], basex=10, basey=10, color="#FF75FF", linestyle="-", marker="v", markerfacecolor="None")

plt.legend()
plt.xlabel(r"$log_{10}h$")
plt.ylabel(r"$log_{10}E_r$")
plt.ylim(1e-4,4e-1)
plt.xlim(1e-2, 3e-1)

plt.savefig("Figure3.png")
