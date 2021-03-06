from asyncore import read
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_theme(style="darkgrid")
# plt.rcParams['text.usetex'] = True
DATA_DIR = r"../data/"
PLOTS_DIR = r"."
FIGURE_FORMAT='pdf'

def savefig(name):
    plt.savefig('{}/{}.{}'.format(PLOTS_DIR,name,FIGURE_FORMAT),bbox_inches='tight')

def read_data(name):
    return pd.read_csv(DATA_DIR+name)

PEAK_PERFORMANCE = 118.4
DGEMM_PERFORMANCE = 82.8702

q3 = read_data("q3.csv")
plt.plot(q3["m"],q3["GFLOPS_mean"],label="Mean GFLOPS")
plt.plot(q3["m"],(q3["m"]==q3["m"])*PEAK_PERFORMANCE,label="Peak Performance")
plt.plot(q3["m"],(q3["m"]==q3["m"])*DGEMM_PERFORMANCE,label="DGEMM reference")
lower_bound = q3["GFLOPS_mean"] - q3["GFLOPS_std"]
upper_bound = q3["GFLOPS_mean"] + q3["GFLOPS_std"]

plt.fill_between(q3["m"], lower_bound, upper_bound, alpha=.3)
plt.xlabel("m")
plt.ylim(bottom=0)
plt.ylabel("GFLOPS")
plt.title(r"GFLOPS performed for m in [1e3,1e4]")
plt.legend(bbox_to_anchor=(1.05, 1))
savefig("q3")
plt.show()

q3 = read_data("q3-50.csv")
plt.plot(q3["m"],q3["GFLOPS_mean"],label="Mean GFLOPS")
plt.plot(q3["m"],(q3["m"]==q3["m"])*PEAK_PERFORMANCE,label="Peak Performance")
plt.plot(q3["m"],(q3["m"]==q3["m"])*DGEMM_PERFORMANCE,label="DGEMM reference")
lower_bound = q3["GFLOPS_mean"] - q3["GFLOPS_std"]
upper_bound = q3["GFLOPS_mean"] + q3["GFLOPS_std"]

plt.fill_between(q3["m"], lower_bound, upper_bound, alpha=.3)
plt.xlabel("m")
plt.ylim(bottom=0)
plt.ylabel("GFLOPS")
plt.title(r"GFLOPS performed for m in [1e3,1e4]")
plt.legend(bbox_to_anchor=(1.05, 1))
savefig("q3-50")
plt.show()

q4 = read_data("q4.csv")
plt.plot(q4["n"],q4["time"][0]/q4["time"],label="Actual Speedup")
plt.plot(q4["n"],q4["n"],label="Ideal Speedup")
plt.xlabel("# number of cores")
plt.ylabel("Speedup")
plt.title(r"Speedup of procedure for different levels of OpenMP parallelization")
plt.legend()
savefig("q4-speedup")
plt.show()

plt.plot(q4["n"],q4["time"][0]/(q4["time"]*q4["n"]),label="Actual Efficiency")
plt.plot(q4["n"],q4["n"]==q4["n"],label="Ideal Efficiency")
plt.xlabel("number of cores")
plt.ylabel("Efficiency")
plt.title(r"Efficiency of procedure for different levels of OpenMP parallelization")
plt.legend()
savefig("q4-efficiency")
plt.show()

