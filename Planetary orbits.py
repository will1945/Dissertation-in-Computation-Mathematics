import matplotlib.pyplot as plt
from math import pi
from math import sin
from math import cos
from math import sqrt


#Initial conditions for numerical scheme
a = GM = 1        #Assumption given
e = 0.5           #Eccentricity (measure of the epllipticity of the orbit)
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
P = 2*pi                             #Period of the orbit
h = P/1000                           #Time step size 
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate
star_position = [0,0]

#Initial conditions for analytical solution
b = sqrt(1-e**2)

xan = []
yan = []
for i in range(0,101):
    E = (2*pi/100)*i
    xan.append(cos(E)+e)
    yan.append(b*sin(E))
    

#Forward Euler
xnum = []
ynum = []
for j in range(0,1000):
    r_n1 = [x + h*vx, y + h*vy]                 #New position vector
    r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))     #New distance from star to planet
    F_r = (-x/r_n0**3,-y/r_n0**3)
    v_n1 = [vx + h*F_r[0],vy + h*F_r[1]]
    x = r_n1[0]
    xnum.append(x)
    y = r_n1[1]
    ynum.append(y)
    vx = v_n1[0]
    vy = v_n1[1]

     
plt.figure("Planetary orbit using Forward Euler scheme")
plt.title("Planetary orbit using Forward Euler scheme")
plt.plot(star_position[0],star_position[1],marker='o', markerfacecolor='yellow', markersize=10)
plt.plot(xnum,ynum,label = "Forward Euler")
plt.plot(xan,yan,label = "Analytical")
plt.legend()
plt.show()



#Modified Euler ---------------------------------------------------------------
#I.C's
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate

xnum_mod = []
ynum_mod = []
for j in range(0,1000):
    r_n1_m = [x + h*vx, y + h*vy]            #New position vector
    d1 = sqrt(((star_position[0]-r_n1_m[0])**2 + (star_position[1]-r_n1_m[1])**2))     #New distance from star to planet
    x = r_n1_m[0]
    y = r_n1_m[1]
    xnum_mod.append(x)
    ynum_mod.append(y)
    F_r1_m = (-x/d1**3,-y/d1**3)
    v_n1_m = [vx + h*F_r1_m[0],vy + h*F_r1_m[1]]
    vx = v_n1_m[0]
    vy = v_n1_m[1]

     
plt.figure("Planetary orbit using Modified Euler scheme")
plt.title("Planetary orbit using Modified Euler scheme")
plt.plot(star_position[0],star_position[1],marker='o', markerfacecolor='yellow', markersize=10)
plt.plot(xnum_mod,ynum_mod,label = "Modified Euler")
plt.plot(xan,yan,label = "Analytical")
plt.legend()
plt.show()
    

    
#Leapfrog ----------------------------------------------------------------------
#I.C's
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate

xleapfrog = []
yleapfrog = []
for j in range(0,1000):
    dr = [x + h/2*vx, y + h/2*vy]            #New position vector
    d = sqrt((star_position[0]-dr[0])**2 + (star_position[1]-dr[1])**2)     #New distance from star to planet
    F_r_leap = [-dr[0]/d**3,-dr[1]/d**3]
    v_1_leap = [vx + h*F_r_leap[0],vy + h*F_r_leap[1]]
    r_1_leap = [dr[0]+h/2*v_1_leap[0],dr[1]+h/2*v_1_leap[1]]
    x = r_1_leap[0]
    y = r_1_leap[1]
    xleapfrog.append(x)
    yleapfrog.append(y)
    vx = v_1_leap[0]
    vy = v_1_leap[1]

plt.figure("Planetary orbit using Leapfrog scheme")
plt.title("Planetary orbit using Leapfrog scheme")
plt.plot(star_position[0],star_position[1],marker='o', markerfacecolor='yellow', markersize=10)
plt.plot(xleapfrog,yleapfrog,label = "Leapfrog")
plt.plot(xan,yan,label = "Analytical")
plt.legend()
plt.show()



#RK4 ------------------------------------------------------------------------
#I.C's
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate  
w = [x,y,vx,vy]     
x_RK4 = []
y_RK4 = []
h = P/250

for j in range(0,250):
    d = sqrt((star_position[0]-w[0])**2 + (star_position[1]-w[1])**2)   
    k1 = [h*w[2],h*w[3],-h*w[0]/d**3,-h*w[1]/d**3]
    
    d = sqrt((star_position[0]-(w[0]+0.5*k1[0]))**2 + (star_position[1]-(w[1]+0.5*k1[1]))**2)  
    k2 = [h*(w[2]+0.5*k1[2]),h*(w[3]+0.5*k1[3]),h*(-(w[0]+0.5*k1[0])/d**3),h*(-(w[1]+0.5*k1[1])/d**3)]
    
    d = sqrt((star_position[0]-(w[0]+0.5*k2[0]))**2 + (star_position[1]-(w[1]+0.5*k2[1]))**2)
    k3 = [h*(w[2]+0.5*k2[2]),h*(w[3]+0.5*k2[3]),h*(-(w[0]+0.5*k2[0])/d**3),h*(-(w[1]+0.5*k2[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+k3[0]))**2 + (star_position[1]-(w[1]+k3[1]))**2)
    k4 = [h*(w[2]+k3[2]),h*(w[3]+k3[3]),h*(-(w[0]+k3[0])/d**3),h*(-(w[1]+k3[1])/d**3)]
   
    w_1 = [w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]),w[2] + (1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]),w[3] + (1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])] 
    w = w_1

    x_RK4.append(w[0])
    y_RK4.append(w[1])
    

plt.figure("Planetary orbit using Runge-Kutta scheme")
plt.title("Planetary orbit using Runge-Kutta scheme")
plt.plot(star_position[0],star_position[1],marker='o', markerfacecolor='yellow', markersize=10)
plt.plot(x_RK4,y_RK4,label = "RK4")
plt.plot(xan,yan,label = "Analytical")
plt.legend()
plt.show()


#Forward Euler energy: 
e = 0.5
h = P/300
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))             
E0 = 0.5*(vx**2+vy**2) - 1/r_n0

time1 = []
frac_E_error1 = []
for j in range(0,300*N):
    r_n1 = [x + h*vx, y + h*vy]              
    r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))     
    F_r = (-x/r_n0**3,-y/r_n0**3)
    v_n1 = [vx + h*F_r[0],vy + h*F_r[1]]
    x = r_n1[0]
    y = r_n1[1]
    vx = v_n1[0]
    vy = v_n1[1]
    E = 0.5*(vx**2+vy**2)-1/r_n0
    frac_E_error1.append(abs((E-E0))/abs(E0))   
    time1.append(j)
 
plt.figure("Energy fractional error")
plt.title("Energy fractional error")
plt.loglog(time1,frac_E_error1,label = "Forward Euler")
plt.xlabel("Time steps")
plt.ylabel("Log((E-E0)/E0)")


#Modified Euler Energy --------------------------------------------------------
e = 0.5
h = P/300
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))             
E0 = 0.5*(vx**2+vy**2) - 1/r_n0

time2 = []
frac_E_error2 = []
for j in range(0,300*N):
    r_n1_m = [x + h*vx, y + h*vy]           
    d1 = sqrt(((star_position[0]-r_n1_m[0])**2 + (star_position[1]-r_n1_m[1])**2))    
    x = r_n1_m[0]
    y = r_n1_m[1]
    F_r1_m = (-x/d1**3,-y/d1**3)
    v_n1_m = [vx + h*F_r1_m[0],vy + h*F_r1_m[1]]
    vx = v_n1_m[0]
    vy = v_n1_m[1]
        
    E = 0.5*(vx**2+vy**2)-1/d1
    frac_E_error2.append(abs((E-E0))/abs(E0))   
    time2.append(j)    
plt.loglog(time2,frac_E_error2,label = "Modified Euler")


#Leapfrog Energy --------------------------------------------------------------
e = 0.5
h = P/300
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))             
E0 = 0.5*(vx**2+vy**2) - 1/r_n0

time3 = []
frac_E_error3 = []
for j in range(0,300*N):
    dr = [x + h/2*vx, y + h/2*vy]           
    d = sqrt((star_position[0]-dr[0])**2 + (star_position[1]-dr[1])**2)    
    F_r_leap = [-dr[0]/d**3,-dr[1]/d**3]
    v_1_leap = [vx + h*F_r_leap[0],vy + h*F_r_leap[1]]
    r_1_leap = [dr[0]+h/2*v_1_leap[0],dr[1]+h/2*v_1_leap[1]]
    x = r_1_leap[0]
    y = r_1_leap[1]
    vx = v_1_leap[0]
    vy = v_1_leap[1]
    E = 0.5*(vx**2+vy**2)-1/(sqrt(x**2+y**2))
    frac_E_error3.append(abs((E-E0))/abs(E0))   
    time3.append(j)     
    
plt.loglog(time3,frac_E_error3,label = "Leapfrog")    
plt.legend()


#RK4 Energy ------------------------------------------------------------------
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate  
w = [x,y,vx,vy]     
E0 = 0.5*(vx**2+vy**2) - 1/sqrt(x**2+y**2)
N = 10**2
h = P/300

frac_E_error4 = []
for j in range(0,300*N):
    d = sqrt((star_position[0]-w[0])**2 + (star_position[1]-w[1])**2) 
    k1 = [h*w[2],h*w[3],-h*w[0]/d**3,-h*w[1]/d**3]

    d = sqrt((star_position[0]-(w[0]+0.5*k1[0]))**2 + (star_position[1]-(w[1]+0.5*k1[1]))**2)  
    k2 = [h*(w[2]+0.5*k1[2]),h*(w[3]+0.5*k1[3]),h*(-(w[0]+0.5*k1[0])/d**3),h*(-(w[1]+0.5*k1[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+0.5*k2[0]))**2 + (star_position[1]-(w[1]+0.5*k2[1]))**2)
    k3 = [h*(w[2]+0.5*k2[2]),h*(w[3]+0.5*k2[3]),h*(-(w[0]+0.5*k2[0])/d**3),h*(-(w[1]+0.5*k2[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+k3[0]))**2 + (star_position[1]-(w[1]+k3[1]))**2)
    k4 = [h*(w[2]+k3[2]),h*(w[3]+k3[3]),h*(-(w[0]+k3[0])/d**3),h*(-(w[1]+k3[1])/d**3)]
 
    w_1 = [w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]),w[2] + (1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]),w[3] + (1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])] 
    w = w_1
    
    E = 0.5*(w[2]**2 +w[3]**2)-1/sqrt((w[0]**2 +w[1]**2))
    error = abs(E-E0)/abs(E0)
    frac_E_error4.append(error)

plt.loglog(frac_E_error4,'r',label = "RK4")    
plt.legend()
plt.show()


#Forward Euler angular momentum ----------------------------------------------
e = 0.5
h = P/300
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
L0 = x*vy - y*vx
time1L = []
frac_E_error1L = []

for j in range(0,300*N):
    r_n1 = [x + h*vx, y + h*vy]              
    r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))     
    F_r = (-x/r_n0**3,-y/r_n0**3)
    v_n1 = [vx + h*F_r[0],vy + h*F_r[1]]
    x = r_n1[0]
    y = r_n1[1]
    vx = v_n1[0]
    vy = v_n1[1]
    L = x*vy - y*vx
    frac_E_error1L.append(abs((L-L0))/abs(L0))   
    time1L.append(j)
 
plt.figure("Angular momentum fractional error")
plt.title("Angular momentum fractional error")
plt.loglog(time1L,frac_E_error1L,label = "Forward Euler")
plt.xlabel("Time")
plt.ylabel("Log((L-L0)/L0)")


#Modified Euler angular momentum 
e = 0.5
h = P/300
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
L0 = x*vy - y*vx

time2L = []
frac_E_error2L = []
for j in range(0,300*N):
    r_n1_m = [x + h*vx, y + h*vy]           
    d1 = sqrt(((star_position[0]-r_n1_m[0])**2 + (star_position[1]-r_n1_m[1])**2))    
    x = r_n1_m[0]
    y = r_n1_m[1]
    F_r1_m = (-x/d1**3,-y/d1**3)
    v_n1_m = [vx + h*F_r1_m[0],vy + h*F_r1_m[1]]
    vx = v_n1_m[0]
    vy = v_n1_m[1]
    L = x*vy - y*vx
    frac_E_error2L.append(abs((L-L0))/abs(L0)) 
    time2L.append(j)

plt.loglog(time2L,frac_E_error2L,label = "Modified Euler")


#Leapfrog angular momentum 
e = 0.5
h = P/300
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))         
L0 = x*vy - y*vx

time3L = []
frac_E_error3L = []
for j in range(0,300*N):
    dr = [x + h/2*vx, y + h/2*vy]           
    d = sqrt((star_position[0]-dr[0])**2 + (star_position[1]-dr[1])**2)    
    F_r_leap = [-dr[0]/d**3,-dr[1]/d**3]
    v_1_leap = [vx + h*F_r_leap[0],vy + h*F_r_leap[1]]
    r_1_leap = [dr[0]+h/2*v_1_leap[0],dr[1]+h/2*v_1_leap[1]]     
    x = r_1_leap[0]
    y = r_1_leap[1]
    vx = v_1_leap[0]
    vy = v_1_leap[1]
    L = x*vy - y*vx
    frac_E_error3L.append(abs((L-L0))/abs(L0))   
    time3L.append(j)     

plt.loglog(time3L,frac_E_error3L,label = "Leapfrog")    


#RK4 angular momentum 
e = 0.5
h = P/300
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
w = [x,y,vx,vy]      
L0 = x*vy - y*vx

time4L = []
frac_E_error4L = []
for j in range(0,300*N):
    d = sqrt((star_position[0]-w[0])**2 + (star_position[1]-w[1])**2) 
    k1 = [h*w[2],h*w[3],-h*w[0]/d**3,-h*w[1]/d**3]

    d = sqrt((star_position[0]-(w[0]+0.5*k1[0]))**2 + (star_position[1]-(w[1]+0.5*k1[1]))**2)  
    k2 = [h*(w[2]+0.5*k1[2]),h*(w[3]+0.5*k1[3]),h*(-(w[0]+0.5*k1[0])/d**3),h*(-(w[1]+0.5*k1[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+0.5*k2[0]))**2 + (star_position[1]-(w[1]+0.5*k2[1]))**2)
    k3 = [h*(w[2]+0.5*k2[2]),h*(w[3]+0.5*k2[3]),h*(-(w[0]+0.5*k2[0])/d**3),h*(-(w[1]+0.5*k2[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+k3[0]))**2 + (star_position[1]-(w[1]+k3[1]))**2)
    k4 = [h*(w[2]+k3[2]),h*(w[3]+k3[3]),h*(-(w[0]+k3[0])/d**3),h*(-(w[1]+k3[1])/d**3)]
 
    w_1 = [w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]),w[2] + (1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]),w[3] + (1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])] 
    w = w_1
    
    L = w[0]*w[3] - w[1]*w[2]
    frac_E_error4L.append(abs((L-L0))/abs(L0))  
    time4L.append(j)   
plt.loglog(time4L,frac_E_error4L,label = "RK4")    
plt.legend()
plt.show()



#Plots for 100 orbits ---------------------------------------------------------
#Forward Euler - all orbits
e = 0.5           #Eccentricity (measure of the epllipticity of the orbit)
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
P = 2*pi                             #Period of the orbit
h = P/300                         #Time step size 
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate
N = 10**2

xnumNP = []
ynumNP = []
for j in range(0,300*N):
    r_n1 = [x + h*vx, y + h*vy]                 #New position vector
    r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))     #New distance from star to planet
    F_r = (-x/r_n0**3,-y/r_n0**3)
    v_n1 = [vx + h*F_r[0],vy + h*F_r[1]]
    x = r_n1[0]
    xnumNP.append(x)
    y = r_n1[1]
    ynumNP.append(y)
    vx = v_n1[0]
    vy = v_n1[1]
     
plt.figure("100 Planetary orbits using Forward Euler scheme")
plt.title("100 Planetary orbits using Forward Euler scheme")
plt.plot(star_position[0],star_position[1],marker='o', markerfacecolor='yellow', markersize=10)
plt.plot(xnumNP,ynumNP,label = "Forward Euler")
plt.plot(xan,yan,label = "Analytical")
plt.legend()
plt.show()


#Modified Euler - all orbits 
e = 0.5           #Eccentricity (measure of the epllipticity of the orbit)
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
P = 2*pi                             #Period of the orbit
h = P/300                         #Time step size 
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate
N = 10**2

xnum_mod = []
ynum_mod = []
for j in range(0,300*N+1):
    r_n1_m = [x + h*vx, y + h*vy]            #New position vector
    d1 = sqrt(((star_position[0]-r_n1_m[0])**2 + (star_position[1]-r_n1_m[1])**2))     #New distance from star to planet
    x = r_n1_m[0]
    y = r_n1_m[1]
    xnum_mod.append(x)
    ynum_mod.append(y)
    F_r1_m = (-x/d1**3,-y/d1**3)
    v_n1_m = [vx + h*F_r1_m[0],vy + h*F_r1_m[1]]
    vx = v_n1_m[0]
    vy = v_n1_m[1]
     
plt.figure("100 Planetary orbits using Modified Euler scheme")
plt.title("100 Planetary orbits using Modified Euler scheme")
plt.plot(star_position[0],star_position[1],marker='o', markerfacecolor='yellow', markersize=10)
plt.plot(xnum_mod,ynum_mod,label = "Modified Euler")
plt.plot(xan,yan,label = "Analytical")
plt.legend()
plt.show()


#Leapfrog - all orbits
e = 0.5           #Eccentricity (measure of the epllipticity of the orbit)
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
P = 2*pi                             #Period of the orbit
h = P/300                         #Time step size 
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate
N = 10**2

xleapfrog = []
yleapfrog = []
for j in range(0,300*N+1):
    dr = [x + h/2*vx, y + h/2*vy]            #New position vector
    d = sqrt((star_position[0]-dr[0])**2 + (star_position[1]-dr[1])**2)     #New distance from star to planet
    F_r_leap = [-dr[0]/d**3,-dr[1]/d**3]
    v_1_leap = [vx + h*F_r_leap[0],vy + h*F_r_leap[1]]
    r_1_leap = [dr[0]+h/2*v_1_leap[0],dr[1]+h/2*v_1_leap[1]]
    x = r_1_leap[0]
    y = r_1_leap[1]
    xleapfrog.append(x)
    yleapfrog.append(y)
    vx = v_1_leap[0]
    vy = v_1_leap[1]
    
plt.figure("100 Planetary orbits using Leapfrog scheme")
plt.title("100 Planetary orbits using Leapfrog scheme")
plt.plot(star_position[0],star_position[1],marker='o', markerfacecolor='yellow', markersize=10)
plt.plot(xleapfrog,yleapfrog,label = "Leapfrog")
plt.plot(xan,yan,label = "Analytical")
plt.legend()
plt.show()


#RK4- all orbits
e = 0.5           #Eccentricity (measure of the epllipticity of the orbit)
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
P = 2*pi                             #Period of the orbit
h = P/75                        #Time step size 
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate
N = 10**2

x_RK4 = []
y_RK4 = []
w = [x,y,vx,vy] 
for j in range(0,75*N+1):
    d = sqrt((star_position[0]-w[0])**2 + (star_position[1]-w[1])**2) 
    k1 = [h*w[2],h*w[3],-h*w[0]/d**3,-h*w[1]/d**3]

    d = sqrt((star_position[0]-(w[0]+0.5*k1[0]))**2 + (star_position[1]-(w[1]+0.5*k1[1]))**2)  
    k2 = [h*(w[2]+0.5*k1[2]),h*(w[3]+0.5*k1[3]),h*(-(w[0]+0.5*k1[0])/d**3),h*(-(w[1]+0.5*k1[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+0.5*k2[0]))**2 + (star_position[1]-(w[1]+0.5*k2[1]))**2)
    k3 = [h*(w[2]+0.5*k2[2]),h*(w[3]+0.5*k2[3]),h*(-(w[0]+0.5*k2[0])/d**3),h*(-(w[1]+0.5*k2[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+k3[0]))**2 + (star_position[1]-(w[1]+k3[1]))**2)
    k4 = [h*(w[2]+k3[2]),h*(w[3]+k3[3]),h*(-(w[0]+k3[0])/d**3),h*(-(w[1]+k3[1])/d**3)]
 
    w_1 = [w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]),w[2] + (1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]),w[3] + (1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])] 
    w = w_1
    
    x_RK4.append(w[0])
    y_RK4.append(w[1])
    

plt.figure("100 planetary orbits using Runge-Kutta scheme")
plt.title("100 planetary orbits using Runge-Kutta scheme")
plt.plot(star_position[0],star_position[1],marker='o', markerfacecolor='yellow', markersize=10)
plt.plot(x_RK4,y_RK4,label = "RK4")
plt.plot(xan,yan,label = "Analytical")
plt.legend()
plt.show()




#Error plots when e = 0.9 and force evaluations = 1000 ------------------------

#Forward Euler energy: 
e = 0.9
h = P/1000
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))             
E0 = 0.5*(vx**2+vy**2) - 1/r_n0
time1 = []
frac_E_error1 = []

for j in range(0,1000*N):
    r_n1 = [x + h*vx, y + h*vy]              
    r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))     
    F_r = (-x/r_n0**3,-y/r_n0**3)
    v_n1 = [vx + h*F_r[0],vy + h*F_r[1]]
    x = r_n1[0]
    y = r_n1[1]
    vx = v_n1[0]
    vy = v_n1[1]
    E = 0.5*(vx**2+vy**2)-1/r_n0
    frac_E_error1.append(abs((E-E0))/abs(E0))   
    time1.append(j)
 
plt.figure("Energy fractional error, using e = 0.9 and force evaluations = 1000")
plt.title("Energy fractional error, using e = 0.9 and force evaluations = 1000")
plt.loglog(time1,frac_E_error1,label = "Forward Euler")
plt.xlabel("Time steps")
plt.ylabel("Log((E-E0)/E0)")


#Modified Euler Energy 
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))             
E0 = 0.5*(vx**2+vy**2) - 1/r_n0

time2 = []
frac_E_error2 = []
for j in range(0,1000*N):
    r_n1_m = [x + h*vx, y + h*vy]           
    d1 = sqrt(((star_position[0]-r_n1_m[0])**2 + (star_position[1]-r_n1_m[1])**2))    
    x = r_n1_m[0]
    y = r_n1_m[1]
    F_r1_m = (-x/d1**3,-y/d1**3)
    v_n1_m = [vx + h*F_r1_m[0],vy + h*F_r1_m[1]]
    vx = v_n1_m[0]
    vy = v_n1_m[1]
    E = 0.5*(vx**2+vy**2)-1/d1
    frac_E_error2.append(abs((E-E0))/abs(E0))   
    time2.append(j)    
    
plt.loglog(time2,frac_E_error2,label = "Modified Euler")


#Leapfrog Energy
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))             
E0 = 0.5*(vx**2+vy**2) - 1/r_n0

time3 = []
frac_E_error3 = []
for j in range(0,1000*N):
    dr = [x + h/2*vx, y + h/2*vy]           
    d = sqrt((star_position[0]-dr[0])**2 + (star_position[1]-dr[1])**2)    
    F_r_leap = [-dr[0]/d**3,-dr[1]/d**3]
    v_1_leap = [vx + h*F_r_leap[0],vy + h*F_r_leap[1]]
    r_1_leap = [dr[0]+h/2*v_1_leap[0],dr[1]+h/2*v_1_leap[1]]
    x = r_1_leap[0]
    y = r_1_leap[1]
    vx = v_1_leap[0]
    vy = v_1_leap[1]
    E = 0.5*(vx**2+vy**2)-1/(sqrt(x**2+y**2))
    frac_E_error3.append(abs((E-E0))/abs(E0))   
    time3.append(j)   
    
plt.loglog(time3,frac_E_error3,label = "Leapfrog")    
plt.legend()


#RK4 Energy
x = (1+e)         #Initial x coordinate position of planet
y = 0             #Initial y coordinate position of planet
vx = 0                               #Velocity in x coordinate
vy = sqrt((1-e)/(1+e))               #Velocity in y coordinate  
w = [x,y,vx,vy]     
E0 = 0.5*(vx**2+vy**2) - 1/sqrt(x**2+y**2)
N = 10**2
h = P/1000

time4 = []
frac_E_error4 = []
for j in range(0,1000*N):
    d = sqrt((star_position[0]-w[0])**2 + (star_position[1]-w[1])**2) 
    k1 = [h*w[2],h*w[3],-h*w[0]/d**3,-h*w[1]/d**3]

    d = sqrt((star_position[0]-(w[0]+0.5*k1[0]))**2 + (star_position[1]-(w[1]+0.5*k1[1]))**2)  
    k2 = [h*(w[2]+0.5*k1[2]),h*(w[3]+0.5*k1[3]),h*(-(w[0]+0.5*k1[0])/d**3),h*(-(w[1]+0.5*k1[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+0.5*k2[0]))**2 + (star_position[1]-(w[1]+0.5*k2[1]))**2)
    k3 = [h*(w[2]+0.5*k2[2]),h*(w[3]+0.5*k2[3]),h*(-(w[0]+0.5*k2[0])/d**3),h*(-(w[1]+0.5*k2[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+k3[0]))**2 + (star_position[1]-(w[1]+k3[1]))**2)
    k4 = [h*(w[2]+k3[2]),h*(w[3]+k3[3]),h*(-(w[0]+k3[0])/d**3),h*(-(w[1]+k3[1])/d**3)]
 
    w_1 = [w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]),w[2] + (1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]),w[3] + (1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])] 
    w = w_1
    
    E = 0.5*(w[2]**2 +w[3]**2)-1/sqrt((w[0]**2 +w[1]**2))
    error = abs(E-E0)/abs(E0)
    frac_E_error4.append(error)
    time4.append(j)
    
plt.loglog(time4,frac_E_error4,label = "RK4")    
plt.legend()
plt.show()


#Forward Euler angular momentum E = 0.9 AND 1000 ORBITS----------------------------------------------
e = 0.9
h = P/1000
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
L0 = x*vy - y*vx
time1L = []
frac_E_error1L = []

for j in range(0,1000*N):
    r_n1 = [x + h*vx, y + h*vy]              
    r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))     
    F_r = (-x/r_n0**3,-y/r_n0**3)
    v_n1 = [vx + h*F_r[0],vy + h*F_r[1]]
    x = r_n1[0]
    y = r_n1[1]
    vx = v_n1[0]
    vy = v_n1[1]
    L = x*vy - y*vx
    frac_E_error1L.append(abs((L-L0))/abs(L0))   
    time1L.append(j)
 
plt.figure("Angular momentum fractional error e = 0.9, using 1000 orbits")
plt.title("Angular momentum fractional error e = 0.9, using 1000 orbits")
plt.loglog(time1L,frac_E_error1L,label = "Forward Euler")
plt.xlabel("Time")
plt.ylabel("Log((L-L0)/L0)")


#Modified Euler angular momentum 
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
L0 = x*vy - y*vx

time2L = []
frac_E_error2L = []
for j in range(0,1000*N):
    r_n1_m = [x + h*vx, y + h*vy]           
    d1 = sqrt(((star_position[0]-r_n1_m[0])**2 + (star_position[1]-r_n1_m[1])**2))    
    x = r_n1_m[0]
    y = r_n1_m[1]
    F_r1_m = (-x/d1**3,-y/d1**3)
    v_n1_m = [vx + h*F_r1_m[0],vy + h*F_r1_m[1]]
    vx = v_n1_m[0]
    vy = v_n1_m[1]
    L = x*vy - y*vx
    frac_E_error2L.append(abs((L-L0))/abs(L0)) 
    time2L.append(j)

plt.loglog(time2L,frac_E_error2L,label = "Modified Euler")


#Leapfrog angular momentum
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))        
L0 = x*vy - y*vx

time3L = []
frac_E_error3L = []
for j in range(0,1000*N):
    dr = [x + h/2*vx, y + h/2*vy]           
    d = sqrt((star_position[0]-dr[0])**2 + (star_position[1]-dr[1])**2)    
    F_r_leap = [-dr[0]/d**3,-dr[1]/d**3]
    v_1_leap = [vx + h*F_r_leap[0],vy + h*F_r_leap[1]]
    r_1_leap = [dr[0]+h/2*v_1_leap[0],dr[1]+h/2*v_1_leap[1]]
    x = r_1_leap[0]
    y = r_1_leap[1]
    vx = v_1_leap[0]
    vy = v_1_leap[1]
    
    L = x*vy - y*vx
    frac_E_error3L.append(abs((L-L0))/abs(L0))   
    time3L.append(j)     

plt.loglog(time3L,frac_E_error3L,label = "Leapfrog")    


#RK4 angular momentum
h = P/1000
N = 10**2
x = (1+e)       
y = 0  
vx = 0                               
vy = sqrt((1-e)/(1+e))   
w = [x,y,vx,vy] 
L0 = x*vy - y*vx

time4L = []
frac_E_error4L = []
for j in range(0,1000*N):
    d = sqrt((star_position[0]-w[0])**2 + (star_position[1]-w[1])**2) 
    k1 = [h*w[2],h*w[3],-h*w[0]/d**3,-h*w[1]/d**3]

    d = sqrt((star_position[0]-(w[0]+0.5*k1[0]))**2 + (star_position[1]-(w[1]+0.5*k1[1]))**2)  
    k2 = [h*(w[2]+0.5*k1[2]),h*(w[3]+0.5*k1[3]),h*(-(w[0]+0.5*k1[0])/d**3),h*(-(w[1]+0.5*k1[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+0.5*k2[0]))**2 + (star_position[1]-(w[1]+0.5*k2[1]))**2)
    k3 = [h*(w[2]+0.5*k2[2]),h*(w[3]+0.5*k2[3]),h*(-(w[0]+0.5*k2[0])/d**3),h*(-(w[1]+0.5*k2[1])/d**3)]

    d = sqrt((star_position[0]-(w[0]+k3[0]))**2 + (star_position[1]-(w[1]+k3[1]))**2)
    k4 = [h*(w[2]+k3[2]),h*(w[3]+k3[3]),h*(-(w[0]+k3[0])/d**3),h*(-(w[1]+k3[1])/d**3)]
 
    w_1 = [w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]),w[2] + (1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]),w[3] + (1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])] 
    w = w_1
    
    L = w[0]*w[3] - w[1]*w[2]
    frac_E_error4L.append(abs((L-L0))/abs(L0))  
    time4L.append(j)   
plt.loglog(time4L,frac_E_error4L,label = "RK4")    
plt.legend()
plt.show()



#Time-reversibility -----------------------------------------------------
#Foward Euler 
e = 0.5          
x = (1+e)        
y = 0            
P = 2*pi                            
h = P/300                        
vx = 0                               
vy = sqrt((1-e)/(1+e))              
N = 10

xnumNP = []
ynumNP = []
for j in range(0,300*N):
    r_n1 = [x + h*vx, y + h*vy]                 #New position vector
    r_n0 = sqrt(((star_position[0]-x)**2 + (star_position[1]-y)**2))     #New distance from star to planet
    F_r = (-x/r_n0**3,-y/r_n0**3)
    v_n1 = [vx + h*F_r[0],vy + h*F_r[1]]
    x = r_n1[0]
    y = r_n1[1]
    vx = v_n1[0]
    vy = v_n1[1]
    
re_x = x
re_y = y
re_vx = -vx
re_vy = -vy

for j in range(0,300*N):
    r_n1 = [re_x + h*re_vx, re_y + h*re_vy]                 #New position vector
    r_n0 = sqrt(((star_position[0]-re_x)**2 + (star_position[1]-re_y)**2))     #New distance from star to planet
    F_r = (-re_x/r_n0**3,-re_y/r_n0**3)
    v_n1 = [re_vx + h*F_r[0],re_vy + h*F_r[1]]
    re_x = r_n1[0]
    re_y = r_n1[1]
    re_vx = v_n1[0]
    re_vy = v_n1[1]
    
final_x = re_x
final_y = re_y
final_vx = -re_vx
final_vy = -re_vy

Euler_Error_x = (1+e) - final_x  
Euler_Error_y = 0 - final_y 
Euler_Error_vx = 0 - final_vx
Euler_Error_vy = sqrt((1-e)/(1+e)) - final_vy
Euler_position_error = sqrt(Euler_Error_x**2 + Euler_Error_y**2)

print("Euler position error=", Euler_position_error)


#Modified Euler
e = 0.5          
x = (1+e)        
y = 0            
P = 2*pi                            
h = P/300                        
vx = 0                               
vy = sqrt((1-e)/(1+e))              
N = 10

xnum_mod = []
ynum_mod = []
for j in range(0,300*N+1):
    r_n1_m = [x + h*vx, y + h*vy]            #New position vector
    d1 = sqrt(((star_position[0]-r_n1_m[0])**2 + (star_position[1]-r_n1_m[1])**2))     #New distance from star to planet
    x = r_n1_m[0]
    y = r_n1_m[1]
    F_r1_m = (-x/d1**3,-y/d1**3)
    v_n1_m = [vx + h*F_r1_m[0],vy + h*F_r1_m[1]]
    vx = v_n1_m[0]
    vy = v_n1_m[1]
    
re_x = x
re_y = y
re_vx = -vx
re_vy = -vy


for j in range(0,300*N+1):
    r_n1_m = [x + h*vx, y + h*vy]            #New position vector
    d1 = sqrt(((star_position[0]-r_n1_m[0])**2 + (star_position[1]-r_n1_m[1])**2))     #New distance from star to planet
    x = r_n1_m[0]
    y = r_n1_m[1]
    F_r1_m = (-x/d1**3,-y/d1**3)
    v_n1_m = [vx + h*F_r1_m[0],vy + h*F_r1_m[1]]
    vx = v_n1_m[0]
    vy = v_n1_m[1]
final_x = re_x
final_y = re_y
final_vx = -re_vx
final_vy = -re_vy

Mod_Euler_Error_x = (1+e) - final_x  
Mod_Euler_Error_y = 0 - final_y 
Mod_Euler_Error_vx = 0 - final_vx
Mod_Euler_Error_vy = sqrt((1-e)/(1+e)) - final_vy
Mod_Euler_position_error = sqrt(Mod_Euler_Error_x**2 + Mod_Euler_Error_y**2)

print("Modified Euler position error=", Mod_Euler_position_error)


#Leapfrog - in the absence of numerical round-off errors, we should arrive back exactly at the original 
e = 0.5           
x = (1+e)        
y = 0            
P = 2*pi                            
h = P/300                        
vx = 0                               
vy = sqrt((1-e)/(1+e))               
N = 10

for j in range(0,300*N+1):
    dr = [x + h/2*vx, y + h/2*vy]            #New position vector
    d = sqrt((star_position[0]-dr[0])**2 + (star_position[1]-dr[1])**2)     #New distance from star to planet
    F_r_leap = [-dr[0]/d**3,-dr[1]/d**3]
    v_1_leap = [vx + h*F_r_leap[0],vy + h*F_r_leap[1]]
    r_1_leap = [dr[0]+h/2*v_1_leap[0],dr[1]+h/2*v_1_leap[1]]
    x = r_1_leap[0]
    y = r_1_leap[1]
    vx = v_1_leap[0]
    vy = v_1_leap[1]
re_x = x
re_y = y
re_vx = -vx
re_vy = -vy

for j in range(0,300*N+1):
    dr = [re_x + h/2*re_vx, re_y + h/2*re_vy]            #New position vector
    d = sqrt((star_position[0]-dr[0])**2 + (star_position[1]-dr[1])**2)     #New distance from star to planet
    F_r_leap = [-dr[0]/d**3,-dr[1]/d**3]
    v_1_leap = [re_vx + h*F_r_leap[0],re_vy + h*F_r_leap[1]]
    r_1_leap = [dr[0]+h/2*v_1_leap[0],dr[1]+h/2*v_1_leap[1]]
    re_x = r_1_leap[0]
    re_y = r_1_leap[1]
    re_vx = v_1_leap[0]
    re_vy = v_1_leap[1]
final_x = re_x
final_y = re_y
final_vx = -re_vx
final_vy = -re_vy

Leap_Error_x = (1+e) - final_x  
Leap_Error_y = 0 - final_y 
Leap_Error_vx = 0 - final_vx
Leap_Error_vy = sqrt((1-e)/(1+e)) - final_vy
Leap_position_error = sqrt(Leap_Error_x**2 + Leap_Error_y**2)

print("Leap position error=",Leap_position_error)

    