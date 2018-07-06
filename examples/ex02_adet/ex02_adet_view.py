import matplotlib.pylab as plt
import pandas as pd

# read and plot kernel
data = pd.read_csv("kernel.csv", sep=";")
plt.plot(data.x,data.phi, label="A5 kernel")
plt.legend()
plt.savefig("kernel.png")
plt.gcf().clear()

# read and plot results
data = pd.read_csv("results.csv", sep=";")
[ plt.plot(data.x, data[eps], label=eps) for eps in data.filter(regex=("eps*"))]
plt.gca().set_ylim([-2, 2])
plt.legend()
plt.savefig("results.png")
plt.show()

"""
# read and plot cfunc
data = pd.read_csv("cfunc.csv", sep=";")
plt.plot(data.r,data.c, label="lee-yang c-func")
plt.legend()
plt.savefig("cfunc.png")
plt.gcf().clear()
"""
