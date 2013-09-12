import numpy as np
import matplotlib.pylab as plt

infoFile = open("../build-project1-Desktop-Debug/info.txt", "r")
n = int(infoFile.readline())
infoFile.close();

xFile    = open("../build-project1-Desktop-Debug/x_data.txt", "r")
numFile  = open("../build-project1-Desktop-Debug/numerical_data.txt", "r")
analFile = open("../build-project1-Desktop-Debug/analytical_data.txt", "r")

xData    = np.zeros(n)
numData  = np.zeros(n)
analData = np.zeros(n)

for i in range(n):
	xData[i]    = float(xFile.readline())
	numData[i]  = float(numFile.readline())
	analData[i] = float(analFile.readline())

plt.plot(xData, analData, xData, numData)
plt.show()

xFile.close()
numFile.close()
analFile.close()
