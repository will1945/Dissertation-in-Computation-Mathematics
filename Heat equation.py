import numpy as np 
import matplotlib.pyplot as plt
from math import pi
from math import sqrt


L = 1       #Length of metal bar
T = 36      #Time steps
N = 5       #Number of points 
N2 = 100
C = 1/3     #Courant number
K = 1       #Thermal diffusivity 
dx = 1/(N-1)                #Change in x 
dt = C*(dx)**2/K            #Change in t
x = np.linspace(0,L,N)      #Space mesh
xa = np.linspace(0,L,N2)      #Space mesh
Temp_I = np.sin(pi*x)  #Inital temperature
Temp_N = np.zeros(N)   #New temperature


#Compute tempurature values 
for j in range(0,T+1):
    plt.figure("Numerical solution vs Analytical solution") 
    plt.clf()
    plt.plot(x,Temp_I,label = "Numerical solution")     
    plt.legend()
    plt.xlabel("Length")
    plt.ylabel("Temperature")
    plt.title("t = 0.75") 
    
    for i in range(1,N-1):
        Temp_N[i] = Temp_I[i] + C*(Temp_I[i+1]-2*Temp_I[i]+Temp_I[i-1])
    Temp_I[:] = Temp_N        #Numberical solutions    
    t = dt*j  
    Temp_A = np.sin(pi*xa)*np.exp((-K*(pi)**2)*t)     #Analytical solutions 
    
    plt.plot(xa,Temp_A,label = "Analytical solution")
    plt.legend()
    plt.show()


def solve_E(C,N,T):
    dx = 1/(N-1)
    K = 1
    dt = C*(dx)**2/K
    x = np.linspace(0,L,N)          
    Temp_I = np.sin(pi*x)        
    Temp_N = np.zeros(N)         
    for j in range(0,T+1):
        t = dt*j  
        Temp_A = np.sin(pi*x)*np.exp((-K*(pi)**2)*t)
    for j in range(0,T):
        for i in range(1,N-1):
            Temp_N[i] = Temp_I[i] + C*(Temp_I[i+1]-2*Temp_I[i]+Temp_I[i-1])
        Temp_I[:] = Temp_N            
    error = (Temp_A - Temp_I)
    E = sqrt(sum(error**2)/N)
    return E


#N vs E error plot:
y = []                   
n = list(range(3,51))

for N in range(3,51):
    T = int((0.5/(1/3*(1/(1-N))**2)))
    y.append(solve_E(1/3,N,T))
plt.figure("Error curve - N vs E")
plt.plot(n,y)
plt.xlabel("Number of points")
plt.ylabel("Error")
plt.title("Plot of E vs N")
plt.show()
plt.figure("Loglog Error curve - N vs E")
plt.loglog(n,y)
plt.xlabel("Number of points")
plt.ylabel("Error")
plt.title("LogLog plot of E vs N")
plt.show()

gradient_NvsE = np.polyfit(np.log(n), np.log(y),1)[0]  


#N vs C error curve:
z = []                    #values of error at each C
xc = []                   #value of Courant number
for c in range(1,49):
    C = 0.01*c
    T = int(8/C)
    print("C=",C,"T=",T,"E=",solve_E(C,5,T))
    z.append(solve_E(C,5,T))
    xc.append(C)
    
    
plt.figure("Error curve - E vs C")
plt.plot(xc,z)
plt.xlabel("Value of Courant number")
plt.ylabel("Error")
plt.title("Plot of E vs C")
plt.show()

plt.figure("LogLog error curve - E vs C")
plt.loglog(xc,z)
plt.xlabel("Value of Courant number")
plt.ylabel("Error")
plt.title("LogLog plot of E vs C")
plt.show()

gradient_NvsC = np.polyfit(np.log(xc), np.log(z),1)[0]  
print("NvsC grad =",gradient_NvsC) 



#Crank Nicolson Scheme -------------------------------------------------------
T = 12
N = 5
I = np.identity(N-2) 
D2 = np.array([[-2/(dx)**2,1/(dx)**2,0/(dx)**2],[1/(dx)**2,-2/(dx)**2,1/(dx)**2],[0/(dx)**2,1/(dx)**2,-2/(dx)**2]])
A = I/dt - (K/2)*D2
b = I/dt + (K/2)*D2

Temp_New = np.zeros(N)
Temp_In = np.sin(pi*x)  #Inital tempurature

for i in range(0,T+1):
    plt.figure("Crank Nicolson Scheme") 
    plt.plot(x,Temp_In,label = "Numberical solution")        
    plt.axis([0,L,0,1])
    plt.xlabel("Length")
    plt.ylabel("Temperature")
    Temp_In = np.array([Temp_In[1],Temp_In[2],Temp_In[3]])
    bu = np.dot(b,Temp_In)
    Temp_New[1:-1] = (np.linalg.solve(A,bu))
    Temp_In = Temp_New


def solve_Crank_E(C,N,T):
    dx = 1/(N-1)
    K = 1
    dt = C*(dx)**2/K
    x = np.linspace(0,L,N)    
    Temp_In = np.sin(pi*x)   
    Temp_New = np.zeros(N) 
    I = np.identity(N-2) 
    D2 = np.array([[-2/(dx)**2,1/(dx)**2,0/(dx)**2],[1/(dx)**2,-2/(dx)**2,1/(dx)**2],[0/(dx)**2,1/(dx)**2,-2/(dx)**2]])
    A = I/dt - (K/2)*D2
    b = I/dt + (K/2)*D2   
    for j in range(0,T+1):
        t = dt*j  
        Temp_A = np.sin(pi*x)*np.exp((-K*(pi)**2)*t)  
    for j in range(0,T+1):
        Temp_In = np.array([Temp_In[1],Temp_In[2],Temp_In[3]])
        bu = np.dot(b,Temp_In)
        Temp_New[1:-1] = (np.linalg.solve(A,bu))
        Temp_In = Temp_New 
    error = (Temp_A - Temp_In)
    E = sqrt(sum(error**2)/N)
    return E


#Crank Nicolson scheme N vs C error graph:
zcn = []                    #values of error at each C
xcn = []                   #value of Courant number
for c in range(39,201):
    C = 0.01*c
    T = int(8/C)
    zcn.append(solve_Crank_E(C,5,T))
    xcn.append(C)
    
plt.figure("Crank Nicholson - E vs C")
plt.plot(xcn,zcn)
plt.xlabel("Value of Courant number")
plt.ylabel("Error")
plt.title("Plot of E vs C")
plt.show()    

plt.figure("Loglog Crank Nicolson - E vs C")
plt.loglog(xcn,zcn)
plt.xlabel("Value of Courant number")
plt.ylabel("Error")
plt.title("LogLog Plot of E vs C")
plt.show()    
gradient_CRANK = np.polyfit(np.log(xcn), np.log(zcn),1)[0]

    

