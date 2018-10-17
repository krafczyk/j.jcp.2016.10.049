import csv
import matplotlib.pyplot as plt
import math


plt.figure(3, figsize=(12,12), dpi=320)

# Load Figure 3 nonconvex data
nonconvex_data = []
first = True
with open("Figure3_nonconvex.out", "r") as nonconvex_datafile:
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
with open("Figure3_convex.out", "r") as convex_datafile:
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

def make_log10(alist):
    return [math.log10(x) for x in alist]

def make_slope2line(x1, x2, y1):
    y2 = 2*(x2-x1)+y1
    return {'x':[x1,x2],'y':[y1,y2]}

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# gamma = 0.25 plots
slopeline_range = [min(make_log10(nonconvex_data_dict[0.25][3.0]['h'])), max(make_log10(nonconvex_data_dict[0.25][3.0]['h']))]

plt.subplot(321)

plt.plot([], color="#FF0000", linestyle="-", marker="x", label="$\\tau=0.505$")
plt.plot([], color="#02FF02", linestyle="-", marker="o", markerfacecolor="None", label="$\\tau=0.6$")
plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")
plt.plot([], color="#FF75FF", linestyle="-", marker="s", markerfacecolor="None", label="$\\tau=3$")
plt.plot([], color="#000000", linestyle="--", label="Slope=2")

plt.plot(make_log10(nonconvex_data_dict[0.25][0.505]['h']), make_log10(nonconvex_data_dict[0.25][0.505]['e']), color="#FF0000", linestyle="-", marker="x")
plt.plot(make_log10(nonconvex_data_dict[0.25][0.6]['h']), make_log10(nonconvex_data_dict[0.25][0.6]['e']), color="#02FF02", linestyle="-", marker="o", markerfacecolor="None")
plt.plot(make_log10(nonconvex_data_dict[0.25][1.0]['h']), make_log10(nonconvex_data_dict[0.25][1.0]['e']), color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.plot(make_log10(nonconvex_data_dict[0.25][2.0]['h']), make_log10(nonconvex_data_dict[0.25][2.0]['e']), color="#000000", linestyle="-", marker="v", markerfacecolor="None")
plt.plot(make_log10(nonconvex_data_dict[0.25][3.0]['h']), make_log10(nonconvex_data_dict[0.25][3.0]['e']), color="#FF75FF", linestyle="-", marker="v", markerfacecolor="None")
slopeline = make_slope2line(slopeline_range[0], slopeline_range[1], -4)
plt.plot(slopeline['x'], slopeline['y'], color="#000000", linestyle="--")

plt.legend(frameon=False)
plt.xlabel(r"$\log_{10}h$")
plt.ylabel(r"$\log_{10}E_r$")
plt.ylim(-4,-0.7)
plt.xlim(-2, -0.9)

plt.subplot(322)

plt.plot([], color="#FF0000", linestyle="-", marker="x", label="$\\tau=0.505$")
plt.plot([], color="#02FF02", linestyle="-", marker="o", markerfacecolor="None", label="$\\tau=0.6$")
plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")
plt.plot([], color="#FF75FF", linestyle="-", marker="s", markerfacecolor="None", label="$\\tau=3$")
plt.plot([], color="#000000", linestyle="--", label="Slope=2")

plt.plot(make_log10(convex_data_dict[0.25][0.505]['h']), make_log10(convex_data_dict[0.25][0.505]['e']), color="#FF0000", linestyle="-", marker="x")
plt.plot(make_log10(convex_data_dict[0.25][0.6]['h']), make_log10(convex_data_dict[0.25][0.6]['e']), color="#02FF02", linestyle="-", marker="o", markerfacecolor="None")
plt.plot(make_log10(convex_data_dict[0.25][1.0]['h']), make_log10(convex_data_dict[0.25][1.0]['e']), color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.plot(make_log10(convex_data_dict[0.25][2.0]['h']), make_log10(convex_data_dict[0.25][2.0]['e']), color="#000000", linestyle="-", marker="v", markerfacecolor="None")
plt.plot(make_log10(convex_data_dict[0.25][3.0]['h']), make_log10(convex_data_dict[0.25][3.0]['e']), color="#FF75FF", linestyle="-", marker="v", markerfacecolor="None")
slopeline = make_slope2line(slopeline_range[0], slopeline_range[1], -4)
plt.plot(slopeline['x'], slopeline['y'], color="#000000", linestyle="--")

plt.legend(frameon=False)
plt.xlabel(r"$\log_{10}h$")
plt.ylabel(r"$\log_{10}E_r$")
plt.ylim(-4,-0.5)
plt.xlim(-2, -0.9)

# gamma = 0.75 plots
slopeline_range = [min(make_log10(nonconvex_data_dict[0.75][2.0]['h'])), max(make_log10(nonconvex_data_dict[0.75][2.0]['h']))]

plt.subplot(323)

plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")
plt.plot([], color="#000000", linestyle="--", label="Slope=2")

plt.plot(make_log10(nonconvex_data_dict[0.75][1.0]['h']), make_log10(nonconvex_data_dict[0.75][1.0]['e']), color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.plot(make_log10(nonconvex_data_dict[0.75][2.0]['h']), make_log10(nonconvex_data_dict[0.75][2.0]['e']), color="#000000", linestyle="-", marker="v", markerfacecolor="None")
slopeline = make_slope2line(slopeline_range[0], slopeline_range[1], -3.8)
plt.plot(slopeline['x'], slopeline['y'], color="#000000", linestyle="--")

plt.legend(frameon=False)
plt.xlabel(r"$\log_{10}h$")
plt.ylabel(r"$\log_{10}E_r$")
plt.ylim(-4,-0.7)
plt.xlim(-2, -0.9)

plt.subplot(324)

plt.plot([], color="#FF0000", linestyle="-", marker="x", label="$\\tau=0.505$")
plt.plot([], color="#02FF02", linestyle="-", marker="o", markerfacecolor="None", label="$\\tau=0.6$")
plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")
plt.plot([], color="#FF75FF", linestyle="-", marker="s", markerfacecolor="None", label="$\\tau=3$")
plt.plot([], color="#000000", linestyle="--", label="Slope=2")

plt.plot(make_log10(convex_data_dict[0.75][0.505]['h']), make_log10(convex_data_dict[0.75][0.505]['e']), color="#FF0000", linestyle="-", marker="x")
plt.plot(make_log10(convex_data_dict[0.75][0.6]['h']), make_log10(convex_data_dict[0.75][0.6]['e']), color="#02FF02", linestyle="-", marker="o", markerfacecolor="None")
plt.plot(make_log10(convex_data_dict[0.75][1.0]['h']), make_log10(convex_data_dict[0.75][1.0]['e']), color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.plot(make_log10(convex_data_dict[0.75][2.0]['h']), make_log10(convex_data_dict[0.75][2.0]['e']), color="#000000", linestyle="-", marker="v", markerfacecolor="None")
plt.plot(make_log10(convex_data_dict[0.75][3.0]['h']), make_log10(convex_data_dict[0.75][3.0]['e']), color="#FF75FF", linestyle="-", marker="v", markerfacecolor="None")
slopeline = make_slope2line(slopeline_range[0], slopeline_range[1], -4.3)
plt.plot(slopeline['x'], slopeline['y'], color="#000000", linestyle="--")

plt.legend(frameon=False)
plt.xlabel(r"$\log_{10}h$")
plt.ylabel(r"$\log_{10}E_r$")
plt.ylim(-4.3,0)
plt.xlim(-2, -0.9)

# gamma = 1.0 plots
slopeline_range = [min(make_log10(nonconvex_data_dict[1.0][1.0]['h'])), max(make_log10(nonconvex_data_dict[1.0][1.0]['h']))]

plt.subplot(325)

plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="--", label="Slope=2")

plt.plot(make_log10(nonconvex_data_dict[1.0][1.0]['h']), make_log10(nonconvex_data_dict[1.0][1.0]['e']), color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
slopeline = make_slope2line(slopeline_range[0], slopeline_range[1], -3.8)
plt.plot(slopeline['x'], slopeline['y'], color="#000000", linestyle="--")

plt.legend(frameon=False)
plt.xlabel(r"$\log_{10}h$")
plt.ylabel(r"$\log_{10}E_r$")
plt.ylim(-4,-0.7)
plt.xlim(-2, -0.95)

plt.subplot(326)

plt.plot([], color="#FF0000", linestyle="-", marker="x", label="$\\tau=0.505$")
plt.plot([], color="#02FF02", linestyle="-", marker="o", markerfacecolor="None", label="$\\tau=0.6$")
plt.plot([], color="#0000FF", linestyle="-", marker="v", markerfacecolor="None", label="$\\tau=1$")
plt.plot([], color="#000000", linestyle="-", marker="d", markerfacecolor="None", label="$\\tau=2$")
plt.plot([], color="#FF75FF", linestyle="-", marker="s", markerfacecolor="None", label="$\\tau=3$")
plt.plot([], color="#000000", linestyle="--", label="Slope=2")

plt.plot(make_log10(convex_data_dict[1.0][0.505]['h']), make_log10(convex_data_dict[1.0][0.505]['e']), color="#FF0000", linestyle="-", marker="x")
plt.plot(make_log10(convex_data_dict[1.0][0.6]['h']), make_log10(convex_data_dict[1.0][0.6]['e']), color="#02FF02", linestyle="-", marker="o", markerfacecolor="None")
plt.plot(make_log10(convex_data_dict[1.0][1.0]['h']), make_log10(convex_data_dict[1.0][1.0]['e']), color="#0000FF", linestyle="-", marker="v", markerfacecolor="None")
plt.plot(make_log10(convex_data_dict[1.0][2.0]['h']), make_log10(convex_data_dict[1.0][2.0]['e']), color="#000000", linestyle="-", marker="v", markerfacecolor="None")
plt.plot(make_log10(convex_data_dict[1.0][3.0]['h']), make_log10(convex_data_dict[1.0][3.0]['e']), color="#FF75FF", linestyle="-", marker="v", markerfacecolor="None")
slopeline = make_slope2line(slopeline_range[0], slopeline_range[1], -4.3)
plt.plot(slopeline['x'], slopeline['y'], color="#000000", linestyle="--")

plt.legend(frameon=False)
plt.xlabel(r"$\log_{10}h$")
plt.ylabel(r"$\log_{10}E_r$")
plt.ylim(-4.3,0)
plt.xlim(-2, -0.9)

plt.subplots_adjust(left=0.07, bottom=0.05, right=0.98, top=0.98)

plt.savefig("Figure3.png")
