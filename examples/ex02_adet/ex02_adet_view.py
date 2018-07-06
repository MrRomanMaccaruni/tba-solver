import matplotlib.pylab as plt
import pandas as pd

# read and plot kernel
data = pd.read_csv("kernel.csv", sep=";")
plt.plot(data.x,data.phi, label="A5 kernel")
plt.grid(linestyle=':')
plt.legend()
plt.savefig("kernel.png")
plt.gcf().clear()

# read and plot results
data = pd.read_csv("results.csv", sep=";")
fig, (ax1, ax2) = plt.subplots(1,2,figsize=(7, 4))
[ax1.plot(data.x, data[eps], label=eps) for eps in data.filter(regex=("eps.*"))]
ax1.set_ylim([-2, 0.1])
ax1.grid(linestyle=':')
ax1.legend()
[ax2.plot(data.x, data[L], label=L) for L in data.filter(regex=("L.*"))]
ax2.set_ylim([-0.1, 2])
ax2.grid(linestyle=':')
ax2.legend()
plt.savefig("results.png")
plt.gcf().clear()

"""
# read and plot cfunc
data = pd.read_csv("cfunc.csv", sep=";")
plt.plot(data.r,data.c, label="lee-yang c-func")
plt.legend()
plt.savefig("cfunc.png")
plt.gcf().clear()
"""
