# Figure for linear relative permeability curves

import matplotlib.pyplot as plt

s = [0, 0.1, 0.9, 1]
kr = [0, 0, 1, 1]

plt.plot(s, kr)
plt.xlabel('$s_p$')
plt.ylabel('$k^r_p$')
plt.savefig('relative_permeability_linear.svg')
plt.savefig('relative_permeability_linear.pdf')
