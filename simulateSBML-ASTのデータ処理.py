
import pandas as pd
import matplotlib.pyplot as plt
import libsbml as lib

res_csv = pd.read_csv('./res.csv')
column_size = len(res_csv)
reader = lib.SBMLReader()

if reader == None:
  print("We failed to read SBML")  

doc = reader.readSBMLFromFile("tyson.xml")
model = doc.getModel()
los = model.getListOfSpecies()
los_size = los.size()
IdVector = []

for i in range(los_size):
  Id = los[i].getId()
  IdVector.append(Id)

fig, ax = plt.subplots(facecolor="w")

time_column = res_csv[res_csv.keys()[0]]
for i in range(1,los_size+1):
  each_molecule_column = res_csv[res_csv.keys()[i]]
  plt.xlabel(res_csv.keys()[0])
  plt.ylabel(res_csv.keys()[i])
  ax.plot(time_column, each_molecule_column, linestyle='solid', marker='.', label = IdVector[i-1],lw=1)
ax.set_xlabel("time /[s]")
ax.set_ylabel("molecule density")
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
plt.show()
