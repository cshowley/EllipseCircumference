import numpy as np
import matplotlib.pyplot as plt

#Trapezoidal Rule
def trap(f,a,b,s): # (lambda x:x*np.exp(x),0,1,number of sections)
    #integral = []
    value = []
    step_size = []
    for n in range(1,s):        
        h = float(b - a) / n
        #print h
        intsum = 0
        intsum += f(a) / 2.0
        for i in range (a,n):
            y = (f(a + i * h))
            intsum += y
            #integral.append(y)
        intsum += f(b) / 2.0
        answer = np.float32(intsum * h)
        value.append(answer)
        step_size.append(h)
            
    real_int = []
    for i in range(a,n):
        j = 1.0
        real_int.append(j)
    #error = [abs(x1-x2) for (x1,x2) in zip(integral,real_int)]
    error = [abs(x1-x2)/x2 for (x1,x2) in zip(value,real_int)]
    #print error
    #print step_size
    plt.figure(1)
    plt.title('Plotting Error in Calculating $xe^{x}$')
    plt.loglog(step_size,error,'.')
    plt.xlabel(r'log$\Delta x$')
    plt.ylabel(r'log$\varepsilon$')
    plt.grid(True)
    #Exact solution to integral is 3/4

def ellipse(s):
    eccentricity = []
    circum = []
    for i in range(0,101):
        e = i * .01
        print e
        eccentricity.append(e)
        #e = np.sqrt(1-B**2/A**2)
        for n in range(11,1000):
            #Calculate perimeter with n steps
            h = (np.pi / 2) / n
            intsum = (h / 3) * np.sqrt(1 - e**2 * (np.sin(0))**2)
            for j in range(1,n,2):
                y = 4 * (h / 3) * np.sqrt(1 - e**2 * (np.sin(j * h))**2)
                intsum += y
            for j in range(2,n-1,2):
                y = 2 * (h / 3) * np.sqrt(1 - e**2 * (np.sin(j * h))**2)
                intsum += y
            intsum += (h / 3) * np.sqrt(1 - e**2 * (np.sin(np.pi / 2))**2)
            integral1 = intsum
            #Calculate perimeter with n-1 steps
            h = (np.pi / 2) / (n - 1)
            intsum = (h / 3) * np.sqrt(1 - e**2 * (np.sin(0))**2)
            for j in range(1,n-1,2):
                y = 4 * (h / 3) * np.sqrt(1 - e**2 * (np.sin(j * h))**2)
                intsum += y
            for j in range(2,n-2,2):
                y = 2 * (h / 3) * np.sqrt(1 - e**2 * (np.sin(j * h))**2)
                intsum += y
            intsum += (h / 3) * np.sqrt(1 - e**2 * (np.sin(np.pi / 2))**2)
            integral2 = intsum
            if abs(integral2 - integral1) < 0.01:
                circum.append(integral1)
                break
    plt.figure(2)
    plt.plot(eccentricity,circum,'.')
    plt.grid(True)
    plt.xlabel('Eccentricity')
    plt.ylabel('Elliptic Integral of the Second Kind')
    plt.title('Variations in Perimeter of an Ellipse as a function of Eccentricity')
    
def derive(s):
    stepSize = []
    errorVals = []
    for p in range(1,6):
        n = 4**p + 1
        h = (2 * np.pi) / n
        stepSize.append(h)
        #Computed derivative
        
        cCos = []
        f = (-3.0/2*np.sin(0) + 2*np.sin(h) - 1.0/2*np.sin(2*h)) / h
        cCos.append(f)
        
        for i in range(1,n):
            f = (-1.0/2 * np.sin((i-1)*h) + 1.0/2*np.sin((i+1)*h)) / h
            cCos.append(f)
        
        f = (3.0/2 * np.sin(2*np.pi) - 2 * np.sin(2*np.pi - h) + 1.0/2 * np.sin(2*np.pi - 2*h)) / h
        cCos.append(f)
    
        #Theoretical derivative
        tCos = []
    
        for i in range(0,n):
            f = np.cos(i * h)
            tCos.append(f)
    
    
        difference = [abs(z2 - z1) for (z1,z2) in zip(cCos,tCos)]
        absolute = [abs(y1) for y1 in tCos]
        error = max(difference) / max(absolute)
        errorVals.append(error)
    plt.figure(3)
    plt.title(r'Finite Difference Method for $\frac{d}{dx}sin(x)$')
    plt.loglog(stepSize,errorVals,'.')
    plt.loglog([(2.0 * np.pi)/11.0],[0.0967561312729],'k*',markersize=10)
    plt.loglog([(2.0 * np.pi)/37.0],[0.00951574575413],'s',markersize=8)
    plt.legend(('Test Points','10% Accuracy','1% Accuracy'), loc='upper left')
    #plt.legend(('Single Precision','Double Precision','Modified Single Precision'), loc='lower left')
    plt.ylabel(r'log $\varepsilon$')
    plt.xlabel('log h')
    
    plt.grid(True)
    
def derive2(s): #To determine the number of steps for accuracy levels
    for n in range(5,s):
        h = (2 * np.pi) / n
        #Computed derivative
        
        cCos = []
        f = (-3.0/2*np.sin(0) + 2*np.sin(h) - 1.0/2*np.sin(2*h)) / h
        cCos.append(f)
        
        for i in range(1,n):
            f = (-1.0/2 * np.sin((i-1)*h) + 1.0/2*np.sin((i+1)*h)) / h
            cCos.append(f)
        
        f = (3.0/2 * np.sin(2*np.pi) - 2 * np.sin(2*np.pi - h) + 1.0/2 * np.sin(2*np.pi - 2*h)) / h
        cCos.append(f)
    
        #Theoretical derivative
        tCos = []
    
        for i in range(0,n):
            f = np.cos(i * h)
            tCos.append(f)
    
    
        difference = [abs(z2 - z1) for (z1,z2) in zip(cCos,tCos)]
        absolute = [abs(y1) for y1 in tCos]
        error = max(difference) / max(absolute)
        #if error < 0.1:
        #    print n
        #    print error
        #    break
        #if error < 0.01:
        #    print n
        #    print error
        #    break
    plt.figure(3)
    plt.loglog([(2.0 * np.pi)/11.0, (2.0 * np.pi)/37.0],[0.0967561312729, 0.00951574575413],'+')