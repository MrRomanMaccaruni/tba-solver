import matplotlib.pylab as plt
import pandas as pd

# read and plot kernel
data = pd.read_csv("kernel.csv", sep=";")
plt.plot(data.x,data.phi, label="lee-yang kernel")
plt.legend()
plt.savefig("kernel.png")
plt.gcf().clear()

# read and plot results
data = pd.read_csv("results.csv", sep=";")
plt.plot(data.x,data.eps, label="lee-yang eps")
plt.plot(data.x,data.L, label="lee-yang L")
plt.gca().set_ylim([-0.1, 1])
plt.legend()
plt.savefig("results.png")
plt.gcf().clear()

# read and plot cfunc
data = pd.read_csv("cfunc.csv", sep=";")
plt.plot(data.r,data.c, label="lee-yang c-func")
plt.legend()
plt.savefig("cfunc.png")
plt.gcf().clear()

