import os
os.chdir("/Users/MRisk/Desktop/SpaceFillingCurves/SpaceFilling")
import BaseFunctions as SFC
import numpy as np
import matplotlib.pyplot as plt

plotA = plt.figure(figsize=(7,7))
plotA = plt.plot(*zip(*SFC.Hilbert_Polygon(5)), color = "black", linewidth = 5)
plotA = plt.axis('off', emit=False)
plt.show()

x1, y1 = zip(*SFC.Hilbert_Polygon(2))
x2, y2 = zip(*SFC.Hilbert_Polygon(3))

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
axes[0].plot(x1, y1, linewidth = 3, color = "black")
axes[1].plot(x2, y2, linewidth = 3, color = "black")
axes[0].grid(b = True)
axes[1].grid(b = True)
axes[0].set_xticks([0, 0.25, 0.5, 0.75, 1])
axes[0].set_yticks([0, 0.25, 0.5, 0.75, 1])
axes[1].set_xticks([0, 0.25, 0.5, 0.75, 1])
axes[1].set_yticks([0, 0.25, 0.5, 0.75, 1])
axes[0].set_xticklabels([0, "", "", "", 1])
axes[0].set_yticklabels([0, "", "", "", 1])
axes[1].set_xticklabels([0, "",  "", "", 1])
axes[1].set_yticklabels([0, "", "", "", 1])
axes[0].set_xlim([0,1])
axes[0].set_ylim([0,1])
axes[1].set_xlim([0,1])
axes[1].set_ylim([0,1])
axes[0].set_xlabel("Generation of Hilbert's Curve, n=2")
axes[1].set_xlabel("Generation of Hilbert's Curve, n=3")
fig.tight_layout()
plt.show()


it1 = SFC.Polya_Dissection(9, np.arctan(4/3))
plotA = plt.figure(figsize=(20,10))
for i in it1:
    plotA = plt.plot(*zip(*i), color = "orange", linewidth = 3)
plotA = plt.plot(*zip(*SFC.Polya_Polygon(9, np.arctan(4/3))), color="blue", linewidth=1)
plotA = plt.axis('off', emit=False)
#plotA = plt.xlabel("Generation of Polya's Space Filling Curve, 11th Iteration")
plt.show()


Polya1 = SFC.Polya_Dissection(1, np.arctan(4/3))
Polya2 = SFC.Polya_Dissection(2, np.arctan(4/3))
Polya3 = SFC.Polya_Dissection(3, np.arctan(4/3))
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(9, 1))
for i in Polya1:
    axes[0].plot(*zip(*i), linewidth = 5, color = "black")
axes[0].plot(*zip(*SFC.Polya_Polygon(1, np.arctan(4/3))), color="red", linewidth=3)
for i in Polya2:
    axes[1].plot(*zip(*i), linewidth = 5, color = "black")
axes[1].plot(*zip(*SFC.Polya_Polygon(2, np.arctan(4/3))), color="red", linewidth=3)
for i in Polya3:
    axes[2].plot(*zip(*i), linewidth = 5, color = "black")
axes[2].plot(*zip(*SFC.Polya_Polygon(3, np.arctan(4/3))), color="red", linewidth=3)
axes[0].set_xticks([])
axes[0].set_yticks([])
axes[1].set_xticks([])
axes[1].set_yticks([])
axes[2].set_xticks([])
axes[2].set_yticks([])
axes[0].set_xticklabels([])
axes[0].set_yticklabels([])
axes[1].set_xticklabels([])
axes[1].set_yticklabels([])
axes[2].set_xticklabels([])
axes[2].set_yticklabels([])
axes[0].set_xlim([-0.25,2.25])
axes[0].set_ylim([-0.25,1.25])
axes[1].set_xlim([-0.25,2.25])
axes[1].set_ylim([-0.25,1.25])
axes[2].set_xlim([-0.25,2.25])
axes[2].set_ylim([-0.25,1.25])
axes[0].set_xlabel("Generation of Polya's Curve, n=1")
axes[1].set_xlabel("Generation of Polya's Curve, n=2")
axes[2].set_xlabel("Generation of Polya's Curve, n=3")
fig.tight_layout()
plt.show()

Polya0 = SFC.Polya_Dissection(0, np.arctan(4/3))
Polya2 = SFC.Polya_Dissection(2, np.arctan(4/3))
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 2))
for i in Polya0:
    axes[0].plot(*zip(*i), linewidth = 2, color = "black")
axes[0].arrow(0.15, 0.05, 1.65, 0, width = 0.007, head_width = 0.035, head_length = 0.05, color='black')
for i in Polya2:
    axes[1].plot(*zip(*i), linewidth = 2, color = "black")
axes[1].arrow(0.12, 0.05, 0.44, 0, width=0.007, head_width=0.035, head_length=0.05, color='black')
axes[1].arrow(0.67, 0.1, 0, 0.7, width=0.0125, head_width=0.05, head_length=0.035, color='black')
axes[1].arrow(0.77, 0.85, 0, -0.7, width=0.0125, head_width=0.05, head_length=0.035, color='black')
axes[1].arrow(0.83, 0.05, 1, 0, width=0.007, head_width=0.035, head_length=0.05, color='black')
axes[0].annotate(r"$\theta$", (0.06, 0.025), fontsize=10)
axes[0].annotate("$(0,0)$", (-0.1, -0.1), fontsize=10)
axes[0].annotate("$(x,y)$", (0.7, 1), fontsize=10)
axes[0].annotate("$(2,0)$", (2, -0.1), fontsize=10)
axes[1].annotate(r"$\theta$", (0.06, 0.025), fontsize=10)
axes[1].annotate("$(0,0)$", (-0.1, -0.1), fontsize=10)
axes[1].annotate("$(x,y)$", (0.7, 1), fontsize=10)
axes[1].annotate("$(2,0)$", (2, -0.1), fontsize=10)
axes[1].annotate("1", (0.27, 0.15), weight='bold', fontsize=10)
axes[1].annotate("2", (0.52, 0.36), weight='bold', fontsize=10)
axes[1].annotate("3", (0.92, 0.59), weight='bold', fontsize=10)
axes[1].annotate("4", (1.23, 0.245), weight='bold', fontsize=10)
axes[0].set_xticks([])
axes[0].set_yticks([])
axes[1].set_xticks([])
axes[1].set_yticks([])
axes[0].set_xticklabels([])
axes[0].set_yticklabels([])
axes[1].set_xticklabels([])
axes[1].set_yticklabels([])
axes[0].set_xlim([-0.25,2.25])
axes[0].set_ylim([-0.25,1.25])
axes[1].set_xlim([-0.25,2.25])
axes[1].set_ylim([-0.25,1.25])
axes[0].set_xlabel("Generation of Polya's Curve, n=0")
axes[1].set_xlabel("Generation of Polya's Curve, n=2")
fig.tight_layout()
plt.show()