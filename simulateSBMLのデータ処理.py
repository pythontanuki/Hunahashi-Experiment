
import pandas as pd
import matplotlib.pyplot as plt

res_csv = pd.read_csv('./res.csv')
time_column = res_csv[res_csv.keys()[0]]
s1_column = res_csv[res_csv.keys()[1]]
s2_column = res_csv[res_csv.keys()[2]]
fig, ax = plt.subplots(facecolor="w")
plt.xlabel(res_csv.keys()[0])
plt.ylabel(res_csv.keys()[1])
ax.plot(time_column, s1_column, linestyle='solid', marker='o', label = "s1")
ax.plot(time_column, s2_column, linestyle='solid', marker='o', label = "s2")
ax.set_xlabel("time /[s]")
ax.set_ylabel("molecule density")
ax.legend()
plt.show()
