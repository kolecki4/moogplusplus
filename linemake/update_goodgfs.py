import os

import numpy as np
import pandas as pd

def key(i):
    return float(i[0:12])

with open("goodgf.old") as f:
    goodgfs = f.readlines()

newData = pd.read_csv("C1-VALD.csv")    
for i in range(len(newData)):
    goodgfs.append("%10.3f%10.1f%10.3f%10.4g%40s\n" % (newData["WL_air(A)"][i], 6.0, newData["Excit(eV)"][i],newData["log gf"][i], "VALD3"))
    
    
newData = pd.read_csv("Mg1-VALD.csv")
for i in range(len(newData)):
    goodgfs.append("%10.3f%10.1f%10.3f%10.4g%40s\n" % (newData["WL_air(A)"][i], 12.0, newData["Excit(eV)"][i],newData["log gf"][i], "VALD3"))
        
newData = pd.read_csv("Al1-VALD.csv")
for i in range(len(newData)):
    goodgfs.append("%10.3f%10.1f%10.3f%10.4g%40s\n" % (newData["WL_air(A)"][i], 13.0, newData["Excit(eV)"][i],newData["log gf"][i], "VALD3"))
       
    
   
goodgfs.sort(key = key)



with open("goodgf", "w") as f:
    f.writelines(goodgfs)
