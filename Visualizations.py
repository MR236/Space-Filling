import os
os.chdir("/Users/MRisk/Desktop/SpaceFillingCurves/SpaceFilling")
import numpy as np
import BaseFunctions as SFC
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

# plotA = plt.figure(figsize=(7,7))
# plotA = plt.plot(*zip(*SFC.Sierpinski_Polygon(7)), color = "black", linewidth = 2)
# #plotA = plt.axis('off', emit=False)
# plt.show()

it1 = SFC.Polya_Dissection(5, [0.75, 1])
print(it1)

plotA = plt.figure(figsize=(7,7))
for i in it1:
    plotA = plt.plot(*zip(*i), color = "black", linewidth = 2)
#plotA = plt.axis('off', emit=False)
plt.show()