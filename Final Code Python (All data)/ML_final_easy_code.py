
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 10:19:38 2021

@author: Shivang Srivastava
"""

# Importing the libraries
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from sklearn import *

# Importing the dataset
plt. clf() 
dataset = pd.read_csv(r"C:\Users\ThinkPad\Videos\AllData.csv");
deltaA = np.array(dataset.iloc[:, 0].values)
deltaA=deltaA.reshape(-1,1)
deltaC= np.array(dataset.iloc[:, 1].values)
deltaC=deltaC.reshape(-1,1)
force= np.array(dataset.iloc[:, 2].values)
force=force.reshape(-1,1)
capdata= np.array(dataset.iloc[:, 3].values)
capdata=capdata.reshape(-1,1)
X=(np.concatenate((deltaA,deltaC,force,capdata),axis=1))


#Create ( first guess) : output (Average)
listavg=[]
j=-1;
for items in deltaA:
    j+=1
    listavg.append(0.5*deltaA[j]+0.5*deltaC[j])
    
capdata_avg=np.asarray(listavg)
capdata_avg = capdata_avg.reshape(-1, 1)
y=capdata_avg

#creating test and training dataset
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 1/3, random_state = 0)

#for more use in future
#import ML model (many others in the previous file and you can contact me for all ML codes)
from sklearn.ensemble import RandomForestRegressor
regressor = RandomForestRegressor(n_estimators = 5, random_state = 0)
regressor.fit(X_train, y_train)

#predict model

j=-1
list1=[]
list2=[]
for items in y_test:
    j+=1;
    X_test1=X_test[j,:]
    X_test1=X_test1.reshape(1,-1)
    list1.append(regressor.predict((X_test1)))
  #  print("Prediction of Test Set Data using Random Forest Regressor: "+str((regressor.predict((X_test1)))))
  #  print("Actual test data: "+str(y_test[j,:]))
  #  list2.append(regressor.predict(y_test[j,:]))



plt.scatter(y_test,np.array(list1))   

x = list(range(2, 23))
yy=[]
x=np.array(x)
for items in x:
    
    y =4.0334*math.exp(0.0808*items)
    yy.append(y)
yy=np.array(yy)
#yy=math.exp(x)
plt.plot(x,yy)
plt.xlabel('Experimental average of ΔC/C0 from 2 axes')
plt.ylabel('Predicted ΔC/C0')
plt.title("Prediction vs experimental average")
plt.text(13,7.5,"y = 4.0334*exp(0.0808*x)")
plt.text(13,6.5,"R² = 0.9188")
#plt.xlim([2,23])
plt.show()
#print("\n \n \n \n \n \n \n \n")
print("\n\nx-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x")
print("\n\nFinal Equation:\n\ny = 4.0334 * exp ( 0.0808 * x )\n\n*y is the desired ΔC/C0 and x is the experimental averge of ΔC/C0 from 2 axes*")
print("\n\nx-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x")

print("\n\nPlease input 4 experimental averages of ΔC/C0 from 2 axes on the next line\n")
lst=[]
for i in range(0, 4):
    ele = int(input())
    lst.append(ele) # adding the element
ls2t=[]
for i in range(0, 4):
    ele=4.0334*math.exp(0.0808*float(lst[i]))
    ls2t.append(ele) # adding the element
print("\n\nx-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x")
  
print("\n\nPredicted Flat ΔC/C0's according to reverse engineered function are:")

for itmes in ls2t:
    print("\n\n"+str(itmes))
#print("\n\n:)\n\n")
print("\n\nx-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x")
import sys
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
