import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as scp

#load data
with open('NIPA 1.1.3.csv') as f:
    reader = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    
    date = next(reader),
    data = {  
            
    "gdp" : next(reader),
    "c" : next(reader),
    "investment" : next(reader),
    "fixed_i" : next(reader),
    "nonresidential" : next(reader),
    "structures" : next(reader),
    "equipment" : next(reader),
    "intellectual" : next(reader),
    "residential" : next(reader),
    
    }
f.close()

#take logs
for key in data:
    data[key] = np.log(data[key])
        
#filter - using matrix form from https://www.stata.com/manuals13/tstsfilterhp.pdf
l = 1600 #filter parameter
T = len(data["gdp"])
M = np.zeros((T,T))
M[0][0] = 1 + l
M[T - 1][T - 1] = 1 + l
M[1][1] = 1 + 5*l
M[T - 2][T - 2] = 1 + 5*l
M[0][1] = -2*l
M[1][0] = -2*l
M[T - 2][T - 1] = -2*l
M[T - 1][T - 2] = -2*l
for i in range(2,T - 2):
    M[i][i] = 1 + 6*l
for i in range(1, T - 2):
    M[i][i + 1] = -4*l
    M[i + 1][i] = -4*l
for i in range(2,T):    
    M[i][i - 2] = l
    M[i - 2][i] = l
fil = np.subtract(np.identity(T), np.linalg.inv(M))

c = {} #cyclical component
for key in data:
    c[key] = np.matmul(fil,data[key])
    
plt.figure(1, figsize=(12,9))  
plt.title("GDP (blue), Consumption (orange), Investment (green)")  
plt.plot(np.transpose(date),c["gdp"])
plt.plot(np.transpose(date),c["c"])
plt.plot(np.transpose(date),c["investment"])

corr = {} #correlation with gdp
sd = {} #standard deviation relative to gdp
for key in c:
    corr[key] = scp.pearsonr(c["gdp"],c[key])[0]
    sd[key] = np.std(c[key]) / np.std(c["gdp"])

plt.figure(2, figsize=(12,9))
plt.title("Correlation with GDP")
plt.bar(corr.keys(), corr.values())

plt.figure(3, figsize=(12,9))
plt.title("SD relative to GDP")
plt.bar(sd.keys(), sd.values())

ld = {} #log differences of original series
sdld = {} #standard deviation of ld relative to that of ld["gdp"]
for key in data:
    ld[key] = np.diff(data[key]) 
    sdld[key] = np.std(ld[key]) / np.std(ld["gdp"])
    
plt.figure(4, figsize=(12,9))
plt.title("SD of log changes relative to GDP")
plt.bar(sdld.keys(), sdld.values())

#datet = (np.transpose(date))
#plt.plot(datet[1:len(datet)],ld["gdp"])

autocorr = {} 
autocorrld = {}
for key in data:
    autocorr[key] = np.corrcoef(c[key][0:len(c[key]) - 1], c[key][1:len(c[key])])[0][1]
    autocorrld[key] = np.corrcoef(ld[key][0:len(ld[key]) - 1], ld[key][1:len(ld[key])])[0][1]

plt.figure(5, figsize=(12,9))
plt.title("1 period autocorrelation - deviation from trend")
plt.bar(autocorr.keys(), autocorr.values())  

plt.figure(6, figsize=(12,9))
plt.title("1 period autocorrelation - log changes")
plt.bar(autocorrld.keys(), autocorrld.values())  
    



