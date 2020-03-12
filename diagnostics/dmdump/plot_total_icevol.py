#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# Import required modules.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# --------------------------------------------------------------------------- #

vvl = np.load('vvl.npy')

# --------------------------------------------------------------------------- #

plt.figure()

plt.rc('axes', linewidth=3)

p = plt.plot(np.arange(5*74, 5*147, 5),
             vvl[:, 1]/1.E11,
             '-b', linewidth=3)
ax = plt.gca()
ax.set_aspect(720/2.0E6)


plt.axis([0, 720, 0, 5.1E17])
plt.xticks(np.arange(0, 730, 90))
plt.yticks(np.arange(0., 2.2E6, 0.5E6))
plt.xlabel('Days', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Total Ice Volume', fontname='arial', fontsize=30, fontweight='bold')

ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(18)
    tick.label1.set_fontname('arial')
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(18)
    tick.label1.set_fontname('arial')
    tick.label1.set_fontweight('bold')

# --------------------------------------------------------------------------- #
