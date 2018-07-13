# linear capillary function

import matplotlib.pyplot as plt
import numpy as np

s = np.array([0, 0.1, 0.9, 1])
pc = np.array([-0.1e5, -0.1e5, 0, 0])

plt.plot(s, pc / 1.e3, '-')
plt.xlabel('$s_1$')
plt.ylabel('$P_c$ (kPa)')
plt.savefig('capillary_linear.svg')
plt.savefig('capillary_linear.pdf')
