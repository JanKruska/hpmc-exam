import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_theme(style="darkgrid")
# plt.rcParams['text.usetex'] = True
PLOTS_DIR = r"."
FIGURE_FORMAT='png'

def savefig(name):
    plt.savefig('{}/{}.{}'.format(PLOTS_DIR,name,FIGURE_FORMAT),bbox_inches='tight')


q3 = pd.read_csv("q3.csv")
plt.plot(q3["m"],q3["GFLOPS"])
plt.xlabel("m")
plt.ylabel("GFLOPS")
plt.title(r"GFLOPS performed for m in [1e3,1e4]")
savefig("q3")
plt.show()

q4 = pd.read_csv("q4.csv")
plt.plot(q4["n"],q4["time"][0]/q4["time"],label="Actual Speedup")
plt.plot(q4["n"],q4["n"],label="Ideal Speedup")
# plt.axline([0,0],[1,1])
plt.xlabel("# number of cores")
plt.ylabel("Speedup")
plt.title(r"Speedup of procedure for different levels of OpenMP parallelization")
plt.legend()
savefig("q4")
plt.show()

