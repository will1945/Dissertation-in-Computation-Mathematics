import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gamma

#Initial conditions:
a = np.pi/6
theta = a                             #angle, which satisifies 0 < a < pi
theta_dot = 0
t_max = 10                            #Final time
N = 1000                              #Number of points/time steps
h = t_max/N                           #dt - rescaled nonodimensional timee


#Analytical solution ---------------------------------------------------------
time_values = []
theta_values_anl = []
for i in range(0,N+1):
    time_values.append(i*h)
    theta_values_anl.append(a*np.cos(i*h))
    
plt.figure("Plot: Non-linear, linear and anlytical comparison ")  
plt.plot(time_values,theta_values_anl,'b',label = "Exact linear")
plt.legend()
plt.xlabel("Time (t)")
plt.ylabel(r"$\Theta$")
plt.title("Plot: Non-linear, linear and anlytical solution comparison")
plt.show()

delta_t_values = []
for i in range(0,N+1):
    delta_t_values.append(i*h)

theta_vec_nonlin = [a]
w = (theta,theta_dot)
#Non-linear solution using RK2
for i in range(0,N):
    k1 = (h*w[1],h*-np.sin(w[0])) 
    k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin(w[0]+0.5*k1[0]))
    w1 = (w[0] + k2[0], w[1] + k2[1]) 
    w = w1
    theta_vec_nonlin.append(w1[0])

plt.plot(delta_t_values,theta_vec_nonlin,'r--', label = "Numerical non-linear")
plt.legend()


#Linear solution using RK2 ---------------------------------------------------
theta_vec_lin = [a]
w = (theta,theta_dot)
for i in range(0,N):
    k1 = (h*w[1],h*-w[0]) 
    k2 = (h*(w[1]+0.5*k1[1]),h*-(w[0]+0.5*k1[0]))
    w1 = (w[0] + k2[0], w[1] + k2[1]) 
    w = w1

    theta_vec_lin.append(w1[0])
    
plt.plot(delta_t_values,theta_vec_lin,'y--',label = "Numerical linear")
plt.legend()
plt.show()


#Error for Runge Kuta linear vs analytic solution --------------------------
def ErrorRK2(N):
    theta = a                             
    theta_dot = 0    
    t_max = 10                           
    h = t_max/N      
    anl_sol = a*np.cos(N*h)
    w = (theta,theta_dot)
    for i in range(0,N):
        k1 = (h*w[1],h*-w[0])
        k2 = (h*(w[1]+0.5*k1[1]),h*-(w[0]+0.5*k1[0]))
        w1 = (w[0] + k2[0], w[1] + k2[1]) 
        w = w1       
    lin_sol = w[0]
    E = lin_sol - anl_sol
    return E


def ErrorRK3(N):
    theta = a                             
    theta_dot = 0    
    t_max = 10                           
    h = t_max/N      
    anl_sol = a*np.cos(N*h)
    w = (theta,theta_dot)
    for i in range(0,N):
        k1 = (h*w[1],h*-w[0])
        k2 = (h*(w[1]+k1[1]/3),h*-(w[0]+k1[0]/3))
        k3 = (h*(w[1]+2*k2[1]/3),h*-(w[0]+2*k2[0]/3))
        
        w1 = (w[0] + 1/4*(k1[0]+3*k3[0]),w[1] + 1/4*(k1[1]+3*k3[1]))
        w = w1       
    lin_sol = w[0]
    E = lin_sol - anl_sol
    return E


def ErrorRK4(N):
    theta = a                           
    theta_dot = 0    
    t_max = 10                           
    h = t_max/N      
    anl_sol = a*np.cos(N*h)
    w = (theta,theta_dot)
    for i in range(0,N):
        k1 = (h*w[1],h*-w[0])
        k2 = (h*(w[1]+0.5*k1[1]),h*-(w[0]+0.5*k1[0]))
        k3 = (h*(w[1]+0.5*k2[1]),h*-(w[0]+0.5*k2[0]))
        k4 = (h*(w[1]+k3[1]),h*-(w[0]+k3[0]))
    
        w1 = (w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))
        w = w1       
    lin_sol = w[0]
    E = anl_sol - lin_sol
    return E


E_values_RK2 = []
E_values_RK3 = []
E_values_RK4 = []
dt_values = []
for i in range(1000,11000,1000):
    E_values_RK2.append(ErrorRK2(i))
    E_values_RK3.append(ErrorRK3(i))
    E_values_RK4.append(ErrorRK4(i))
    dt_values.append(10/i)


plt.figure("Loglog plot of Errors Anl vs Linear")    
plt.loglog(dt_values,E_values_RK2,label = "RK2")
plt.loglog(dt_values,E_values_RK3,label = "RK3")
plt.loglog(dt_values,E_values_RK4,label = "RK4")
plt.legend()
plt.xlabel("log h")
plt.ylabel("log E")
plt.title("Loglog plot: E vs h ($\Delta$t) for Runge Kuta schemes")
plt.show()

gradient_RK2 = np.polyfit(np.log(dt_values), np.log(E_values_RK2),1)[0]
print("RK2 =",gradient_RK2)

gradient_RK3 = np.polyfit(np.log(dt_values), np.log(E_values_RK3),1)[0]
print("RK3 =",gradient_RK3)

gradient_RK4 = np.polyfit(np.log(dt_values), np.log(E_values_RK4),1)[0]
print("RK4 =",gradient_RK4)



#Measuring the period --------------------------------------------------------

T_anl = gamma(0.25)**2/np.sqrt(np.pi)   

#Non-Linear 
def Calc_T_RK2_LinearInt(N,a):
    theta = a                             
    theta_dot = 0
    t_max = 10                            
    h = t_max/N                          
    w = (theta,theta_dot)
    thetalist = []
    for i in range(0,N):
        k1 = (h*w[1],h*-np.sin(w[0])) 
        k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin(w[0]+0.5*k1[0]))
        w1 = (w[0] + k2[0], w[1] + k2[1]) 
        w = w1
        thetalist.append(w[0])
        if w[0] < 0:
            x = ((i)*h,(i+1)*h)
            y = (thetalist[-2],thetalist[-1])
            p = np.polyfit(x,y,1)
            linear_root = np.roots(p)
            T = linear_root*4
            break
    return T
  

def Calc_T_RK3_LinearInt(N,a):
    theta = a                             
    theta_dot = 0
    t_max = 10                            
    h = t_max/N                          
    w = (theta,theta_dot)
    thetalist = []
    for i in range(0,N):
        k1 = (h*w[1],h*-np.sin(w[0]))
        k2 = (h*(w[1]+k1[1]/3),h*-np.sin((w[0]+k1[0]/3)))
        k3 = (h*(w[1]+2*k2[1]/3),h*-np.sin((w[0]+2*k2[0]/3)))
        w1 = (w[0] + 1/4*(k1[0]+3*k3[0]),w[1] + 1/4*(k1[1]+3*k3[1]))
        w = w1       
        thetalist.append(w[0])
        if w[0] < 0:
            x = ((i)*h,(i+1)*h)
            y = (thetalist[-2],thetalist[-1])
            p = np.polyfit(x,y,1)
            linear_root = np.roots(p)
            for k in linear_root:
                if (k > 0 and k*4 < t_max):
                    T = k*4
            break
    return T


def Calc_T_RK4_LinearInt(N,a):
    theta = a                            
    theta_dot = 0
    t_max = 10                            
    h = t_max/N                          
    w = (theta,theta_dot)
    thetalist = []
    for i in range(0,N):
        k1 = (h*w[1],h*-np.sin(w[0]))  
        k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin((w[0]+0.5*k1[0])))
        k3 = (h*(w[1]+0.5*k2[1]),h*-np.sin((w[0]+0.5*k2[0])))
        k4 = (h*(w[1]+k3[1]),h*-np.sin((w[0]+k3[0])))
        w1 = (w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))
        w = w1       
        thetalist.append(w[0])
        if w[0] < 0:
            x1 = (i+1)*h
            y1 = w[0]
            x0 = i*h
            y0 = thetalist[-2]
            m = (y1-y0)/(x1-x0)
            T = (x0 - (y0/m))*4
            break
    return T


def Calc_T_RK2_QuadInt(N,a):  
    t_max = 10                           
    h = t_max/N                           
    w = (a,0)
    thetalist = []
    thetalistnew = []
    for i in range(0,N):
        k1 = (h*w[1],h*-np.sin(w[0])) 
        k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin(w[0]+0.5*k1[0]))
        w1 = (w[0] + k2[0], w[1] + k2[1]) 
        w = w1
        thetalist.append(w[0])
        if w[0] < 0:
            n = i
            w = (a,0)
            for j in range(0,n+3):
                k1 = (h*w[1],h*-np.sin(w[0])) 
                k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin(w[0]+0.5*k1[0]))
                w1 = (w[0] + k2[0], w[1] + k2[1]) 
                w = w1
                thetalistnew.append(w[0])
            x = [(n-1)*h,(n)*h,(n+1)*h,(n+2)*h]
            y = [thetalistnew[n-2],thetalistnew[n-1],thetalistnew[n],thetalistnew[n+1]]
            if (abs(y[0]) > abs(y[3])):
                x.remove(x[0])
                y.remove(y[0])
            else:
                x.remove(x[3])
                y.remove(y[3])
            p = np.polyfit(x,y,2)
            quadratic_roots = np.roots(p)
            for k in quadratic_roots:
                if (k > 0 and k*4 < t_max):
                    T = k*4
            break
    return T


def Calc_T_RK2_CubicInt(N,a):
    theta = a                            
    theta_dot = 0
    t_max = 10                          
    h = t_max/N                          
    w = (theta,theta_dot)
    thetalist = []
    thetalistnew = []
    for i in range(0,N):
        k1 = (h*w[1],h*-np.sin(w[0])) 
        k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin(w[0]+0.5*k1[0]))
        w1 = (w[0] + k2[0], w[1] + k2[1]) 
        w = w1
        thetalist.append(w[0])
        if w[0] < 0:
            n = i
            w = (a,0)
            for j in range(0,n+3):
                k1 = (h*w[1],h*-np.sin(w[0])) 
                k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin(w[0]+0.5*k1[0]))
                w1 = (w[0] + k2[0], w[1] + k2[1]) 
                w = w1
                thetalistnew.append(w[0])
            x = [(n-1)*h,(n)*h,(n+1)*h,(n+2)*h]
            y = [thetalistnew[n-2],thetalistnew[n-1],thetalistnew[n],thetalistnew[n+1]]
            p = np.polyfit(x,y,3)
            cubic_roots = np.roots(p)
            for k in cubic_roots:
                if (k > 0 and k*4 < t_max):
                    T = k*4
            break
    return T

    
def Calc_T_RK3_QuadInt(N,a):  
    t_max = 10                      
    h = t_max/N                           
    w = (a,0)
    thetalist = []
    thetalistnew = []
    for i in range(0,N):
        k1 = (h*w[1],h*-np.sin(w[0]))
        k2 = (h*(w[1]+k1[1]/3),h*-np.sin((w[0]+k1[0]/3)))
        k3 = (h*(w[1]+2*k2[1]/3),h*-np.sin((w[0]+2*k2[0]/3)))
        w1 = (w[0] + 1/4*(k1[0]+3*k3[0]),w[1] + 1/4*(k1[1]+3*k3[1]))
        w = w1   
        thetalist.append(w[0])
        if w[0] < 0:
            n = i
            w = (a,0)
            for j in range(0,n+3):
                k1 = (h*w[1],h*-np.sin(w[0]))
                k2 = (h*(w[1]+k1[1]/3),h*-np.sin((w[0]+k1[0]/3)))
                k3 = (h*(w[1]+2*k2[1]/3),h*-np.sin((w[0]+2*k2[0]/3)))
                w1 = (w[0] + 1/4*(k1[0]+3*k3[0]),w[1] + 1/4*(k1[1]+3*k3[1])) 
                w = w1
                thetalistnew.append(w[0])
            x = [(n-1)*h,(n)*h,(n+1)*h,(n+2)*h]
            y = [thetalistnew[n-2],thetalistnew[n-1],thetalistnew[n],thetalistnew[n+1]]
            if (abs(y[0]) > abs(y[3])):
                x.remove(x[0])
                y.remove(y[0])
            else:
                x.remove(x[3])
                y.remove(y[3])
            p = np.polyfit(x,y,2)
            quadratic_roots = np.roots(p)
            for k in quadratic_roots:
                if (k > 0 and k*4 < t_max):
                    T = k*4
            break
    return T


def Calc_T_RK3_CubicInt(N,a):  
    t_max = 10                           
    h = t_max/N                           
    w = (a,0)
    thetalist = []
    thetalistnew = []
    for i in range(0,N):
        k1 = (h*w[1],h*-np.sin(w[0]))
        k2 = (h*(w[1]+k1[1]/3),h*-np.sin((w[0]+k1[0]/3)))
        k3 = (h*(w[1]+2*k2[1]/3),h*-np.sin((w[0]+2*k2[0]/3)))
        w1 = (w[0] + 1/4*(k1[0]+3*k3[0]),w[1] + 1/4*(k1[1]+3*k3[1]))
        w = w1   
        thetalist.append(w[0])
        if w[0] < 0:
            n = i
            w = (a,0)
            for j in range(0,n+2):
                k1 = (h*w[1],h*-np.sin(w[0]))
                k2 = (h*(w[1]+k1[1]/3),h*-np.sin((w[0]+k1[0]/3)))
                k3 = (h*(w[1]+2*k2[1]/3),h*-np.sin((w[0]+2*k2[0]/3)))
                w1 = (w[0] + 1/4*(k1[0]+3*k3[0]),w[1] + 1/4*(k1[1]+3*k3[1]))
                w = w1   
                thetalistnew.append(w[0])
            x = [(n-1)*h,(n)*h,(n+1)*h,(n+2)*h]
            y = [thetalistnew[n-2],thetalistnew[n-1],thetalistnew[n],thetalistnew[n+1]]
            p = np.polyfit(x,y,3)
            cubic_roots = np.roots(p)
            for k in cubic_roots:
                if (k > 0 and k*4 < t_max):
                    T = k*4
            break
    return T


def Calc_T_RK4_QuadInt(N,a):  
    t_max = 10                           
    h = t_max/N                           
    w = (a,0)
    thetalist = []
    thetalistnew = []
    for i in range(0,N):
        k1 = (h*w[1],h*-np.sin(w[0]))  
        k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin((w[0]+0.5*k1[0])))
        k3 = (h*(w[1]+0.5*k2[1]),h*-np.sin((w[0]+0.5*k2[0])))
        k4 = (h*(w[1]+k3[1]),h*-np.sin((w[0]+k3[0])))
        w1 = (w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))
        w = w1 
        thetalist.append(w[0])
        if w[0] < 0:
            n = i
            w = (a,0)
            for j in range(0,n+3):
                k1 = (h*w[1],h*-np.sin(w[0]))  
                k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin((w[0]+0.5*k1[0])))
                k3 = (h*(w[1]+0.5*k2[1]),h*-np.sin((w[0]+0.5*k2[0])))
                k4 = (h*(w[1]+k3[1]),h*-np.sin((w[0]+k3[0])))
                w1 = (w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))
                w = w1 
                thetalistnew.append(w[0])
            x = [(n-1)*h,(n)*h,(n+1)*h,(n+2)*h]
            y = [thetalistnew[n-2],thetalistnew[n-1],thetalistnew[n],thetalistnew[n+1]]
            if (abs(y[0]) > abs(y[3])):
                x.remove(x[0])
                y.remove(y[0])
            else:
                x.remove(x[3])
                y.remove(y[3])
            p = np.polyfit(x,y,2)
            quadratic_roots = np.roots(p)
            for k in quadratic_roots:
                if (k > 0 and k*4 < t_max):
                    T = k*4
            break
    return T


def Calc_T_RK4_CubicInt(N,a):  
    t_max = 10                           
    h = t_max/N                           
    w = (a,0)
    thetalist = []
    thetalistnew = []
    for i in range(0,N):
        k1 = (h*w[1],h*-np.sin(w[0]))  
        k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin((w[0]+0.5*k1[0])))
        k3 = (h*(w[1]+0.5*k2[1]),h*-np.sin((w[0]+0.5*k2[0])))
        k4 = (h*(w[1]+k3[1]),h*-np.sin((w[0]+k3[0])))
        w1 = (w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))
        w = w1 
        thetalist.append(w[0])
        if w[0] < 0:
            n = i
            w = (a,0)
            for j in range(0,n+3):
                k1 = (h*w[1],h*-np.sin(w[0]))  
                k2 = (h*(w[1]+0.5*k1[1]),h*-np.sin((w[0]+0.5*k1[0])))
                k3 = (h*(w[1]+0.5*k2[1]),h*-np.sin((w[0]+0.5*k2[0])))
                k4 = (h*(w[1]+k3[1]),h*-np.sin((w[0]+k3[0])))
                w1 = (w[0] + (1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),w[1] + (1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))
                w = w1 
                thetalistnew.append(w[0])
            x = [(n-1)*h,(n)*h,(n+1)*h,(n+2)*h]
            y = [thetalistnew[n-2],thetalistnew[n-1],thetalistnew[n],thetalistnew[n+1]]
            p = np.polyfit(x,y,3)
            cubic_roots = np.roots(p)
            for k in cubic_roots:
                if (k > 0 and k*4 < t_max):
                    T = k*4
            break
    return T


dt_valuesT = []
T_values_RK2 = []
T_values_RK3 = []
T_values_RK4 = []
T_values_RK2_quad = []
T_values_RK2_cubic = []
T_values_RK3_quad = []
T_values_RK3_cubic = []
T_values_RK4_quad = []
T_values_RK4_cubic = []
for i in range(100,1100,100):
    T_values_RK2.append(abs(Calc_T_RK2_LinearInt(i,np.pi/2)-T_anl))
    T_values_RK3.append(abs(Calc_T_RK3_LinearInt(i,np.pi/2)-T_anl))
    T_values_RK4.append(abs(Calc_T_RK4_LinearInt(i,np.pi/2)-T_anl))
    T_values_RK2_quad.append(abs(Calc_T_RK2_QuadInt(i,np.pi/2)-T_anl))
    T_values_RK2_cubic.append(abs(Calc_T_RK2_CubicInt(i,np.pi/2)-T_anl))
    T_values_RK3_quad.append(abs(Calc_T_RK3_QuadInt(i,np.pi/2)-T_anl))
    T_values_RK3_cubic.append(abs(Calc_T_RK3_CubicInt(i,np.pi/2)-T_anl))
    T_values_RK4_quad.append(abs(Calc_T_RK4_QuadInt(i,np.pi/2)-T_anl))
    T_values_RK4_cubic.append(abs(Calc_T_RK4_CubicInt(i,np.pi/2)-T_anl))
    dt_valuesT.append(10/i)
 

gradient_RK2_T_LIN = np.polyfit(np.log(dt_valuesT), np.log(T_values_RK2),1)[0]
print("RK2 Linear interpolation rate of convergence =",gradient_RK2_T_LIN)

gradient_RK2_T_QUADRATIC = np.polyfit(np.log(dt_valuesT), np.log(T_values_RK2_quad),1)[0]
print("RK2 QUADRATIC rate of convergence =",gradient_RK2_T_QUADRATIC)

gradient_RK2_T_CUBIC = np.polyfit(np.log(dt_valuesT), np.log(T_values_RK2_cubic),1)[0]
print("RK2 CUBIC rate of convergence =",gradient_RK2_T_CUBIC)

gradient_RK3_T_LIN = np.polyfit(np.log(dt_valuesT), np.log(T_values_RK3),1)[0]
print("RK3 Linear interpolation rate of convergence =",gradient_RK3_T_LIN)

gradient_RK3_T_QUAD = np.polyfit(np.log(dt_valuesT), np.log(T_values_RK3_quad),1)[0]
print("RK3 QUADRATIC rate of convergence =",gradient_RK3_T_QUAD)

gradient_RK3_T_CUBIC = np.polyfit(np.log(dt_valuesT), np.log(T_values_RK3_cubic),1)[0]
print("RK3 CUBIC rate of convergence =",gradient_RK3_T_CUBIC)

gradient_RK4_T_LIN = np.polyfit(np.log(dt_valuesT), np.log(T_values_RK4),1)[0]
print("RK4 Linear interpolation rate of convergence =",gradient_RK4_T_LIN)

gradient_RK4_T_QUADRATIC = np.polyfit(np.log(dt_valuesT), np.log(T_values_RK4_quad),1)[0]
print("RK4 QUADRATIC rate of convergence =",gradient_RK4_T_QUADRATIC)

gradient_RK4_T_CUBIC = np.polyfit(np.log(dt_valuesT), np.log(T_values_RK4_cubic),1)[0]
print("RK4 CUBIC rate of convergence =",gradient_RK4_T_CUBIC)


plt.figure("Loglog plot: Rate of convergence RK2")    
plt.title("RK2 rate of convergence")
plt.loglog(dt_valuesT,T_values_RK2,label ="Linear interpolation" )
plt.loglog(dt_valuesT,T_values_RK2_quad,label ="Quadratic interpolation" )
plt.loglog(dt_valuesT,T_values_RK2_cubic,label = "Cubic interpolation")
plt.legend()
plt.xlabel("log h")
plt.ylabel("log $\mid$ Tnum - Tanl $\mid$")
plt.show()


plt.figure("Loglog plot: Rate of convergence RK3")    
plt.title("RK3 rate of convergence")
plt.plot(dt_valuesT,T_values_RK3,label = "Linear interpolation")
plt.plot(dt_valuesT,T_values_RK3_quad,label = "Quadratic interpolation")
plt.plot(dt_valuesT,T_values_RK3_cubic,label = "Cubic interpolation")
plt.xlabel("log h")
plt.ylabel("log $\mid$ Tnum - Tanl $\mid$")
plt.legend()
plt.show()


plt.figure("Loglog plot: Rate of convergence RK4")    
plt.title("RK4 rate of convergence")
plt.plot(dt_valuesT,T_values_RK4,label = "Linear interpolation")
plt.plot(dt_valuesT,T_values_RK4_quad,label = "Quadratic interpolation")
plt.plot(dt_valuesT,T_values_RK4_cubic,label = "Cubic interpolation")
plt.xlabel("log h")
plt.ylabel("log $\mid$ Tnum - Tanl $\mid$")
plt.legend()
plt.show()



#NON-LINEAR -----------------------------------------------------------------------

a_values = [np.pi/6,np.pi/3,np.pi/2,3*np.pi/5,4*np.pi/5]
T_values_vs_a = []
for x in a_values:
    T_values_vs_a.append(abs(Calc_T_RK2_LinearInt(1000,x)))
    
T_values_vs_bigrange_a = []
Tminus2pi_values_vs_bigrange_a = []
avalues = []
for i in range(500,10,-1):
    a = np.pi*(i/1000)
    avalues.append(a)
    T_values_vs_bigrange_a.append(abs(Calc_T_RK2_LinearInt(1000,a)))
    Tminus2pi_values_vs_bigrange_a.append((abs(Calc_T_RK2_LinearInt(1000,a))-2*np.pi))


plt.figure("T vs a")    
plt.title("T vs a")
plt.plot(avalues,T_values_vs_bigrange_a)
plt.xlabel("a")
plt.ylabel("T")
plt.legend()
plt.show() 


plt.figure("T minus 2pi vs a")    
plt.title("T-2$\pi$ vs a")
plt.plot(avalues,Tminus2pi_values_vs_bigrange_a)
plt.xlabel("a ")
plt.ylabel("T-2$\pi$")
plt.legend()
plt.show() 


plt.figure("Loglog T vs a")    
plt.title("T vs a")
plt.loglog(avalues,T_values_vs_bigrange_a)
plt.xlabel("log(a)")
plt.ylabel("log(T)")
plt.legend()
plt.show()  


plt.figure("Loglog T minus 2pi vs a")    
plt.title("T-2$\pi$  vs a")
plt.loglog(avalues,Tminus2pi_values_vs_bigrange_a)
plt.xlabel("log(a)")
plt.ylabel("log(T-2$\pi$)")
plt.legend()
plt.show()  


gradient_T_vs_a= np.polyfit(np.log(avalues), np.log(T_values_vs_bigrange_a),1)[0]
print("RK2 LINEAR T VS a =",gradient_T_vs_a)    

gradient_Tminus2PI_vs_a= np.polyfit(np.log(avalues), np.log(Tminus2pi_values_vs_bigrange_a),1)[0]
print("RK2 LINEAR T-2PI VS a =",gradient_Tminus2PI_vs_a)    