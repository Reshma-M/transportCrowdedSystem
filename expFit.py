import numpy as np
import csv
import itertools
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def read_file(filename):
	data = np.genfromtxt (filename, delimiter="\n")
	#print data
	return data

def read_csvfile(filename):
	with open(filename,'r') as rf:
		reader = csv.reader(rf)
		data = []
		for row in reader:
			data.append([float(i) for i in row])
		#print(np.size(data))
		#data = list(itertools.chain.from_iterable(list(reader)))
		#data = [float(i) for i in data]
		return np.array(data)

def expFit(x, a, b):
    I0 = 0.125; # I0 = 400/8000=0.05, 800/8000=0.1, 1000/8000=0.125
    return I0+(a*(np.exp(-x/b)))

def calLambda(x,inputData):
    popt, pcov = curve_fit(expFit, x, inputData)
    return popt,pcov

n = 20 # sites added
dx=n*8/1000
x1 = np.arange(dx,(dx*10)+dx,dx)
x2 = np.arange(dx,(dx*10)+dx,dx)

# Load data
axData120 = read_csvfile('./densitySummed-20sites-5mins.csv')
axData220 = read_csvfile('./densitySummed-20sites-10mins.csv')


print(np.shape(axData120))
print(np.shape(axData220))

#row-wise mean axis =1
# Extract data near the proximal and distal sites and average over N simulation runs to smooth out the data
distalMean120 = np.mean(np.flip(axData120[:,198-10:198]),axis=0)# data flipped as distal accumulation data is in reverse order of distance from cut site
proximalMean120 = np.mean(axData120[:,203:203+10],axis=0)
distalMean220 = np.mean(np.flip(axData220[:,198-10:198]),axis=0)#988-30:988
proximalMean220 = np.mean(axData220[:,203:203+10],axis=0)#1012:1012+30




allFit = np.zeros((4,6),dtype=float)

print(np.shape(proximalMean120))
print(np.shape(x1))
fit1 = calLambda(x1,proximalMean120)
fit2 = calLambda(x2,distalMean120)
fit3 = calLambda(x1,proximalMean220)
fit4 = calLambda(x2,distalMean220)

allFit[0,:2],allFit[0,2:4],allFit[0,4:] = fit1[0],fit1[1][0],fit1[1][1]
allFit[1,:2],allFit[1,2:4],allFit[1,4:] = fit2[0],fit2[1][0],fit2[1][1]
allFit[2,:2],allFit[2,2:4],allFit[2,4:] = fit3[0],fit3[1][0],fit3[1][1]
allFit[3,:2],allFit[3,2:4],allFit[3,4:] = fit4[0],fit4[1][0],fit4[1][1]

#np.savetxt('./expFit-20sites-5-10mins-I0.csv', allFit, fmt='%10.5f', delimiter=',')

fig,a =  plt.subplots(2,2,figsize=(12, 12))
a[0][0].plot(x1,proximalMean120,'r--',label='prox 5min')
a[0][0].plot(x1, expFit(x1, *fit1[0]), 'c-',label='fit: a=%5.3f, b=%5.3f' % tuple(fit1[0]))
a[0][0].set_title('Proximal 5 min')
a[0][0].legend()
a[0][0].set_xlabel('Distance from cut site (um)')
a[0][0].set_ylabel('Vesicle density')
a[0][1].plot(x1,distalMean120,'r--',label='dist 5min')
a[0][1].plot(x1, expFit(x2, *fit2[0]), 'c-',label='fit: a=%5.3f, b=%5.3f' % tuple(fit2[0]))
a[0][1].set_title('Distal 5 min')
a[0][1].legend()
a[0][1].set_xlabel('Distance from cut site (um)')
a[0][1].set_ylabel('Vesicle density')
a[1][0].plot(x1,proximalMean220,'r--',label='prox 10min')
a[1][0].plot(x1, expFit(x1, *fit3[0]), 'c-',label='fit: a=%5.3f, b=%5.3f' % tuple(fit3[0]))
a[1][0].set_title('Proximal 10 min')
a[1][0].legend()
a[1][0].set_xlabel('Distance from cut site (um)')
a[1][0].set_ylabel('Vesicle density')
a[1][1].plot(x1,distalMean220,'r--',label='dist 10min')
a[1][1].plot(x1, expFit(x2, *fit4[0]), 'c-',label='fit: a=%5.3f, b=%5.3f' % tuple(fit4[0]))
a[1][1].set_title('Distal 10 min')
a[1][1].legend()
a[1][1].set_xlabel('Distance from cut site (um)')
a[1][1].set_ylabel('Vesicle density')

#plt.show()
fig.savefig('./axotomyFit-20sites-5-10min.jpeg')
print(x1)
#fig1, ax1 = plt.subplots()
#ax1.plot(np.mean(axData120[:,238:262],axis=0),label='1min20sec')
#ax1.plot(np.mean(axData220[:,238:262],axis=0),label='2min20sec')
#ax1.plot(np.mean(axData320[:,238:262],axis=0),label='3min20sec')
#ax1.legend()
#fig1.savefig('./../../../cluster/apr2020/axotomy/complexSystem/axotomy-16sites-ves1000.jpeg')
