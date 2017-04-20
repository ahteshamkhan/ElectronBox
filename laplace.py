import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import numpy
import math


#CONSTANT = 10
ELECTRON = -20 #Potential value measured 1 unit from electron
PROTON = 20
K = 1 #Electrostatic constant

class Particle:
    def __init__(self, x0, y0, vx0, vy0, ax0, ay0, mass, charge):
        self.mass = mass
        self.charge = charge

        self.x0, self.y0 = x0, y0
        self.vx0, self.vy0 = vx0, vy0
        self.ax0, self.ay0 = ax0, ay0

        self.r = numpy.array([x0, y0])
        self.v = numpy.array([vx0, vy0])
        self.a = numpy.array([ax0, ay0])

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


class Box:
    """Represents a 2D box whose walls are held at constant potential."""

    def __init__(self, L, W):
        self.length = L
        self.width = W
        self.box = [] #Box matrix
        self.createBox() #Set up box matrix
        self.particles = [] #List of particles in box
        self.particlePositions = []

        #Initialize box boundary potential values
        self.top = 0
        self.bottom = 0
        self.left = 0
        self.right = 0

        "Initialize plot"
        xList = numpy.arange(0, self.width, 1)
        yList = numpy.arange(0, self.length, 1)
        self.X, self.Y = numpy.meshgrid(xList, yList)
        self.Z = numpy.array(self.box)
        fig = plt.figure(figsize=plt.figaspect(0.5))
        self.ax = fig.add_subplot(1, 2, 1, projection='3d')
        
    def setPotentials(self, bottom, top, left, right):
        """Sets the constant potential values of the box."""
        self.top = top
        self.bottom = bottom
        self.left = left
        self.right = right

        # Set the potential values of the four walls
        #print(len(self.box))
        for i in range(self.length):
            self.box[i][0] = self.left
            self.box[i][self.width-1] = self.right
        for i in range(self.width):
            self.box[0][i] = self.top
            self.box[self.length-1][i] = self.bottom

        
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
        i = y
        j = x
        electron = Particle(j, i, 0, 0, 0, 0, 9.11*10**-31, -1.6*10**-19)
        self.particles.append(electron)
        self.particlePositions.append((i, j))
        self.box[i][j] = ELECTRON

    def addProton(self, x, y):
        """Adds an electron to the box matrix."""
        i = y
        j = x
        proton = Particle(j, i, 0, 0, 0, 0, 1.627*10**-27, 1.6*10**-19)
        self.particles.append(proton)
        self.particlePositions.append((i, j))
        self.box[i][j] = PROTON

    def setParticlePositions(self):
        """Sets the particle positions in the box."""
        particlePositions = []
        for particle in self.particles:
            i, j = particle.getR()[1], particle.getR()[0]
            particlePositions.append((int(i), int(j)))
        self.particlePositions = particlePositions

    def relaxation(self, tolerance):
        """Uses the relaxation method to calculate the potential everywhere
           inside the box."""
        #print("Particle positions", self.particlePositions)
        cycle = 0 #Number of full cylces used
        while True:
            #self.surfacePlot2()
            maxChange = 0 #Maximum amount a cell in the box changed values between cycles
            for i in range(len(self.box)):
                cycle += 1
                for j in range(len(self.box[i])):
                    if i != 0 and i != len(self.box)-1 and j != 0 and j != len(self.box[0])-1 and (i, j) not in self.particlePositions:
                        old = self.box[i][j]
                        self.box[i][j] = (self.box[i][j+1] + self.box[i][j-1] + self.box[i+1][j] + self.box[i-1][j])/4
                        if abs(self.box[i][j] - old) > maxChange:
                            maxChange = abs(self.box[i][j] - old)
            if maxChange < tolerance:
                #print("cycle", cycle)
                return

    def electricTensor(self, i, j):
        """Caclulates the electric field at each point."""
        eX = self.box[i][j] - self.box[i][j+1]
        eY = self.box[i][j] - self.box[i+1][j]
        t11 = eX**2 - eY**2
        t12 = eX*eY
        t21 = t12
        t22 = -t11
        return numpy.array([[t11, t12],
                            [t21, t22]])
        

    def maxwellForce(self,k):
        """Uses Maxwell-Stress Tensor to calculate the force on thhe k'th particle."""
        j = int(round(self.particles[k].getR()[0]))
        i = int(round(self.particles[k].getR()[1]))
        #print("pos", i, j)
        force = [0, 0]
        temp = [0, 0]
        i = i - 2
        j = j - 2
        for m in range(4):
            #print(i,j)
            ds = numpy.array([0, 1])
            force += numpy.inner(self.electricTensor(i, j), ds)
            temp += numpy.inner(self.electricTensor(i, j), ds)
            j += 1
        #print(temp)
        temp = [0, 0]
        for m in range(4):
            #print(i,j)
            ds = numpy.array([1, 0])
            force += numpy.inner(self.electricTensor(i, j), ds)
            temp += numpy.inner(self.electricTensor(i, j), ds)
            i += 1
        #print(temp)
        temp = [0, 0]
        for m in range(4):
            #print(i,j)
            ds = numpy.array([0, -1])    
            force += numpy.inner(self.electricTensor(i, j), ds)
            temp += numpy.inner(self.electricTensor(i, j), ds)
            j -= 1
        #print(temp)
        temp = [0, 0]
        for m in range(4):
            #print(i,j)
            ds = numpy.array([-1, 0])
            force += numpy.inner(self.electricTensor(i, j), ds)
            temp += numpy.inner(self.electricTensor(i, j), ds)
            i -= 1
        #print(temp)  
        #print(force)
        return force


    def surfacePlot(self):
        """Rutrns a surface plot of the box."""
        xList = numpy.arange(0, self.width, 1)
        yList = numpy.arange(0, self.length, 1)
        X, Y = numpy.meshgrid(xList, yList)
        Z = numpy.array(self.box)
        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        return ax.plot_wireframe(X, Y, Z, cmap=cm.coolwarm, rstride = 5)

    def surfacePlot2(self):
        """Rutrns a surface plot of the box."""
        
        Z = numpy.array(self.box)
        self.ax.clear()
        self.ax.plot_wireframe(self.X, self.Y, Z, cmap=cm.coolwarm, rstride = 1, cstride = 1)
        plt.pause(.001)

        
    def printBox(self):
        """Prints the box matrix."""
        #print('hello')
        for row in self.box:
            string = ''
            for col in row:
                string += str(col) + ' '
            print(string)

    def eulerStepParticle(self, tau):
        """Positionally advances the particles for one time increment."""
    
        for i in range(len(self.particles)):
            
            self.relaxation(.0000000000000001)
            
            
            accelx = self.maxwellForce(i)[0]/(self.particles[i].getM())
            accely = self.maxwellForce(i)[1]/(self.particles[i].getM())
        
            vx = self.particles[i].getV()[0]
            vy = self.particles[i].getV()[1]
            
            x = self.particles[i].getR()[0]
            y = self.particles[i].getR()[1]
            
            newvx = vx + accelx * tau
            newvy = vy + accely * tau
            
            newx = x + newvx * tau
            newy = y + newvy * tau
            
            self.particles[i].setR([newx, newy])
            self.particles[i].setV([newvx, newvy])
            
            self.setParticlePositions() 
            
            
            
           
        return self.particles


Y = 40
X = 40

box = Box(Y, X)
box.addElectron(18, 18)
#box.addElectron(17, 20)
#box.addElectron(15, 15)
#box.addElectron(20, 20)
#box.addProton(20,20)
#box.addProton(30, 30)

#box.addElectron((int(51/4)+1, int(51/2)+1))

box.setPotentials(0, 1, 0, 0) #Top, bottom, left right 
box.setParticlePositions()
box.relaxation(.0000000000000000001)
#box.relaxation(.001)
#box.surfacePlot2()
#box.relaxation(.00000000000001) #For 31X31
#box.relaxation(.000000000000001) For 51x51

#surface = box.surfacePlot()
print(box.maxwellForce(0))
#print(box.maxwellForce(1))

xpos1 = []
ypos1 = []
xpos2 = []
ypos2 = []
tau = 0.0000000000000005
xpos3 = []
ypos3 = []
xpos4 = []
ypos4 = []

plt.ion()

while(True):
    p = box.eulerStepParticle(tau)
    xpos1.append(p[0].getR()[0])
    ypos1.append(p[0].getR()[1])
#
#    xpos2.append(p[1].getR()[0])
#    ypos2.append(p[1].getR()[1])
#    
#    xpos3.append(p[2].getR()[0])
#    ypos3.append(p[2].getR()[1])
#    xpos4.append(p[3].getR()[0])
#    ypos4.append(p[3].getR()[1])
#    
#    plt.plot(xpos1, ypos1,label =  "electron 1", color = 'red')
#    plt.plot(xpos2, ypos2, label = "electron 2", color = 'orange')
#    plt.plot(xpos3, ypos3, label = "proton 1", color = 'blue')
#    plt.plot(xpos4, ypos4, label = "proton 2", color = 'black')
    box.surfacePlot2()
    plt.pause(0.3)


#surface = box.surfacePlot()
#box.surfacePlot2()
#plt.plot(xpos1, ypos1, color = 'red')
#plt.plot(xpos2, ypos2,  color = 'orange')
#plt.show()
