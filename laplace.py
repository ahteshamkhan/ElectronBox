import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import numpy as np
import math


class Electron:
    def __init__(self, x0, y0, vx0, vy0, ax0, ay0):
        self.mass = 9.10938*10**-31
        self.charge = -1.6*10**-19
            
        self.x0, self.y0 = x0, y0
        self.vx0, self.vy0 = vx0, vy0
        self.ax0, self.ay0 = ax0, ay0
             
        self.r = np.array([x0, y0])
        self.v = np.array([vx0, vy0])
        self.a = np.array([ax0, ay0])
        

    def getR(self):
        return self.r
    
    def getV(self):
        return self.v
    
    def getA(self):
        return self.a
    
    def getM(self):
        return self.mass
    
    def setR(self, r):
        self.r = r
        
    def setV(self, v):
        self.v = v
        
    def setA(self, a):
        self.a = a

class Proton:
    def __init__(self, x0, y0, vx0, vy0):
        self.mass = 1.6726*10**-27
        self.charge = 1.6*10**-19
            
        self.x0, self.y0 = x0, y0
        self.vx0, self.vy0 = vx0, vy0
             
        self.r = np.array([x0, y0])
        self.v = np.array([vx0, vy0])
        
        self.onforce = 0

    def getR(self):
        return self.r
    
    def getV(self):
        return self.v
    
    def getM(self):
        return self.mass
    
    def setR(self, r):
        self.r = r
        
    def setV(self, v):
        self.v = v

    def getonforce(self):
        return self.onforce
    
    def setonforce(self, onforce):
        self.onforce = onforce

#CONSTANT = 10
ELECTRON = -20 #Potential value measured 1 unit from electron
PROTON = 20
#

class Box:
    """Represents a 2D box whose walls are held at constant potential."""

    def __init__(self, L, W):
        self.length = L
        self.width = W
        self.box = [] #Box matrix
        self.createBox() #Set up box matrix
        self.electrons = [] #List of electrons (positions) in box
        self.protons = []

        #Initialize box boundary potential values
        self.top = 0
        self.bottom = 0
        self.left = 0
        self.right = 0

    def setPotentials(self, bottom, top, left, right):
        """Sets the constant potential values of the box."""
        self.top = top
        self.bottom = bottom
        self.left = left
        self.right = right

        # Set the potential values of the four walls
        print("box length", len(self.box))
        for i in range(self.length):
            self.box[i][0] = self.left
            self.box[i][self.width-1] = self.right
        for i in range(self.width):
            self.box[0][i] = self.top
            self.box[self.length-1][i] = self.bottom

        # Set the potential values of each electron
        for electron in self.electrons:
            self.box[electron.getR()[0]][electron.getR()[1]] = ELECTRON
        for proton in self.protons:
            self.box[proton.getR()[0]][proton.getR()[1]] = PROTON
        
    def createBox(self):
        """Creates the box matrix (2-dimensional list)."""
        box = []
        for i in range(self.length):
            row = []
            for j in range(self.width):
                row.append(0)
            box.append(row)
        self.box = box

    def addElectron(self, x, y):
        """Adds an electron to the box matrix."""
        electron = Electron(x, y, 0, 0, 0, 0)
        self.electrons.append(electron)

    def addProton(self, x, y):
        """Adds an electron to the box matrix."""
        proton = Proton(x, y, 0, 0)
        self.protons.append(proton)

    def surfaceCharge(self):
        """Calculates the total charge residing on the inner walls of the box."""

        # Start with potential of 0
        for i in range(1, self.length-1):
            self.box[1][i] = 0
            self.box[self.width-2][i] += 0
        for i in range(1, self.width-1):
            self.box[i][1] += 0
            self.box[i][self.length-2] += 0
        for electron in self.electrons:
            self.box[i][j] = ELECTRON

        # Calculate new potential
        for electron in self.electrons:
            for i in range(1, self.length-1):
                self.box[1][i] += CONSTANT/math.sqrt((1 - electron.getR()[0])**2 + (i - electron.getR()[1])**2)
                self.box[self.width-2][i] += CONSTANT/math.sqrt(((self.width-2) - electron.getR()[0])**2 + (i - electron.getR()[1])**2)
            for i in range(1, self.width-1):
                self.box[i][1] += CONSTANT/math.sqrt((i - electron.getR()[0])**2 + (1 - electron.getR()[1])**2)
                self.box[i][self.length-2] += CONSTANT/math.sqrt((i - electron.getR()[0])**2 + ((self.length-2) - electron.getR()[1])**2)

        # Calculate surface charge
        for i in range(1, self.length-1):
            self.box[0][i] = self.boundary - self.box[1][i]
            self.box[self.width-1][i] = self.boundary - self.box[self.width-2][i]
        for i in range(1, self.width-1):
            self.box[i][0] = self.boundary - self.box[i][1]
            self.box[i][self.length-1] = self.boundary - self.box[i][self.length-2]

    def relaxation(self, tolerance):
        """Uses the relaxation method to calculate the potential everywhere
           inside the box."""
        #print(self.electrons)
        cycle = 0 #Number of full cylces used
        gate = True
        while True:
            maxChange = 0 #Maximum amount a cell in the box changed values between cycles
            for i in range(len(self.box)):
                cycle += 1
                for j in range(len(self.box[i])):
                    if i != 0 and i != len(self.box)-1 and j != 0 and j != len(self.box[0])-1:
                        for electron in range(len(self.electrons)):
                            if i == self.electrons[electron].getR()[0] or j == self.electrons[electron].getR()[1] :
                                gate = False
                        for proton in range(len(self.protons)):
                            if i == self.protons[proton].getR()[0] or j == self.protons[proton].getR()[1] :
                                gate = False
                        if gate == True:
                            old = self.box[i][j]
                            self.box[i][j] = (self.box[i][j+1] + self.box[i][j-1] + self.box[i+1][j] + self.box[i-1][j])/4
                            if abs(self.box[i][j] - old) > maxChange:
                                maxChange = abs(self.box[i][j] - old)
            if maxChange < tolerance:
                print("cycle", cycle)
                return

    def forceE(self, k):
        """Calculates the electric force on the kth electron."""
        i, j = int(self.electrons[k].getR()[0]), int(self.electrons[k].getR()[1])
        return ((self.box[i+1][j]-self.box[i-1][j])/2, (self.box[i][j+1]-self.box[i][j-1])/2)

    def forceP(self, k):
        """Calculates the electric force on the kth proton."""
        i, j = int(self.protons[k].getR()[0]), int(self.protons[k].getR()[1])
        return ((self.box[i+1][j]-self.box[i-1][j])/2, (self.box[i][j+1]-self.box[i][j-1])/2)

    def surfacePlot(self):
        """Rutrns a surface plot of the box."""
        xList = numpy.arange(0, self.width, 1)
        yList = numpy.arange(0, self.length, 1)
        X, Y = numpy.meshgrid(xList, yList)
        Z = numpy.array(self.box)
        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        return ax.plot_wireframe(X, Y, Z, cmap=cm.coolwarm, rstride = 5)

    def printBox(self):
        """Prints the box matrix."""
        #print('hello')
        for row in self.box:
            string = ''
            for col in row:
                string += str(col) + ' '
            print(string)


    def forcesonE(self):
        
        for electron in range(len(self.electrons)):
            thisElectron = self.electrons[electron]
            r = thisElectron.getR()
            v = thisElectron.getV()
            a = thisElectron.getA()
            
            for i in range(1, len(self.electrons)):
                otherElectron = self.electrons[(electron + i)%len(self.electrons)]
                r2 = otherElectron.getR()
                const = 9*10**9 * -1.6*10**-19
                const = const/((np.linalg.norm(r - r2))**3)
                vector = r - r2
                a = a + const*vector
                print(a)
                
            
            for j in range(0, len(self.protons)):
                otherProton = self.protons[(electron + j)%len(self.protons)]
                r2 = otherProton.getR()
                const = 9*10**9 * 1.6*10**-19
                const = const/((np.linalg.norm(r - r2))**3)
                vector = r - r2
                a = a + const*vector
                
            #thisElectron.setA(a)
        
                
    def eulerStepElectron(self, tau):
        """Positionally advances the electrons for one time increment."""
        
        for i in range(len(self.electrons)):
            
            accelx = self.forceE(i)[0]/(9.10938*10**-31)
            accely = self.forceE(i)[1]/(9.10938*10**-31)
            #print(box.electrons[i].getA())
        
            vx = self.electrons[i].getV()[0]
            vy = self.electrons[i].getV()[1]
            
            x = self.electrons[i].getR()[0]
            y = self.electrons[i].getR()[1]
            
            newvx = vx + accelx * tau
            newvy = vy + accely * tau
            
            newx = x + newvx * tau
            newy = y + newvy * tau
            
            self.electrons[i].setR([newx, newy])
            self.electrons[i].setV([newvx, newvy])
            
        return self.electrons
    
    def eulerStepProton(self, tau):
        
        for i in range(len(self.protons)):
            
            accelx = self.forceP(i)[0]/(1.6726*10**-27)
            accely = self.forceP(i)[1]/(1.6726*10**-27)
        
            vx = self.protons[i].getV()[0]
            vy = self.protons[i].getV()[1]
            
            x = self.protons[i].getR()[0]
            y = self.protons[i].getR()[1]
            
            newvx = vx + accelx * tau
            newvy = vy + accely * tau
            
            newx = x + newvx * tau
            newy = y + newvy * tau
            
            self.protons[i].setR([newx, newy])
            self.protons[i].setV([newvx, newvy])
            
        return self.protons

Y = 31
X = 31

box = Box(Y, X)
box.addElectron(int(X/2), int(Y/2))
box.addElectron(16, 20)
box.addProton(25, 18)
box.addProton(2, 10)
#box.addElectron((int(X/2)+1, int(Y/2)))
#box.addElectron((int(51/4)+1, int(51/2)+1))

box.setPotentials(0.1, 0.1, 0.1, 0.1) #Top, bottom, left right 
#box.relaxation(.001) 
box.relaxation(.0000000000000000001) #For 31X31
#box.relaxation(.000000000000001) For 51x51


#for k in range(len(box.electrons)):
#    print(box.forceE(k))
#for k in range(len(box.protons)):
#    print(box.forceP(k))
box.forcesonE()
print(box.electrons[0].getA())
#print(box.forcesonE()[1].getA())

#xpos1 = []
#ypos1 = []
#xpos2 = []
#ypos2 = []
#etstep = 0.000000000000001
#ptstep = 0.0000000000001
#xpos3 = []
#ypos3 = []
#for i in range(0, 30):
#    el = box.eulerStepElectron(etstep)
#    xpos1.append(el[0].getR()[0])
#    ypos1.append(el[0].getR()[1])
#    xpos2.append(el[1].getR()[0])
#    ypos2.append(el[1].getR()[1])
#    
#    pr = box.eulerStepProton(ptstep)
#    xpos3.append(pr[0].getR()[0])
#    ypos3.append(pr[0].getR()[1])
#    
#    
#    
#
#plt.plot(xpos1, ypos1)
#plt.plot(xpos2, ypos2)
#plt.plot(xpos3, ypos3)
surface = box.surfacePlot()
plt.show()
