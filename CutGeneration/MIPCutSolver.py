#!/usr/bin/env python
# coding: utf-8

# In[25]:


import pandas as pd
from gurobipy import *
import os
import openpyxl


# In[ ]:


readPathNoCuts = "./jim_instances/"
readPathCuts = "./CutGenerationLpFiles/"

data = []
for f in os.listdir(readPathNoCuts):
    d = {}
    instanceName = f.split(".")[0]
    d["Instance"] = instanceName
    input1 = (instanceName + ".mps")
    input2 = (instanceName + ".lp")
    m = read(readPathNoCuts + input1)
    m.setParam("Heuristics", 0)
    m.setParam("Cuts", 3)
    m.setParam("NodeLimit", 1)
    m.setParam("TimeLimi", 3600)
    for var in m.getVars():
        var.vType = GRB.INTEGER
    m.update()
    m.optimize()
    d["G+C bound"] = m.ObjBound
    d["G+C obj"] = m.ObjVal
    d["G+C gap"] = str(m.MIPGap) + "%"
    
    m = read(readPathCuts + input2)
    m.setParam("Heuristics", 0)
    m.setParam("Cuts", 3)
    m.setParam("NodeLimit", 1)
    m.setParam("TimeLimit", 3600)
    m.optimize()
    d["G+C+AC bound"] = m.ObjBound
    d["G+C+AC obj"] = m.ObjVal
    d["G+C+AC gap"] = str(m.MIPGap) + "%"
    
    m = read(readPathCuts + input2)
    m.setParam("Heuristics", 0)
    m.setParam("Cuts", 0)
    m.setParam("NodeLimit", 1)
    m.setParam("TimeLimit", 3600)
    m.optimize()
    d["G+AC bound"] = m.ObjBound
    d["G+AC obj"] = m.ObjVal
    d["G+AC gap"] = str(m.MIPGap) + "%"
    


data.append(d)
    
outPut = pd.DataFrame(data)
outPut.to_excel("bounds_objs_and_gaps.xlsx")


# In[ ]:




