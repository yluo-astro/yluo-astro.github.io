import numpy as np
from matplotlib.animation import FuncAnimation 
import matplotlib.pyplot as plt

mh = 1.67262171e-24  #g
me = 9.10938215e-28  #g
pi = 3.14159265358979323846
ev2erg = 1.60217653e-12
hplanck = 6.6260693e-27
kboltz = 1.3806504e-16
clight = 2.99792458e10
StBz = 5.6704e-5
arad = 4*StBz/clight  #7.5657e-15
GravConst = 6.67428e-8 #cm**3 * g**-1 * s**-2
SolarMass=1.9891e33
SolarLuminosity=3.9e33
year = 3.1536e7
pc = 3.0856776e18
au = 1.495978707e13

def set_units(M,L,T):
    mass_unit = M*SolarMass
    length_unit = L*au
    time_unit = T*year
    G_unit = length_unit**3 * mass_unit**-1 * time_unit**-2
    velocity_unit = length_unit/time_unit
    return mass_unit,length_unit,time_unit,G_unit

def get_ic(r0_au,theta0,e):
    r0 = r0_au*au/length_unit
    theta0 = theta0
    eccentricity = e


    ..........
    ..........
    period = .......
    
    return r0, theta0, period, ptheta, E_tot, mass1, mass2, mass_reduced, k_potential
    

def func(r, v, theta):
    r_dot = ..
    v_dot = .......
    theta_dot = ........
    return np.array([r_dot,v_dot,theta_dot])
    
    

def Euler(dt,r,v,theta):
    r_dot,v_dot,theta_dot = func(r,v,theta)
    v_n1 = v + v_dot *dt
    r_n1 = r + r_dot *dt
    theta_n1 = theta + theta_dot*dt
    return r_n1, v_n1, theta_n1
    

def rk4(dt,r,v,theta):    
    k1 = ..
    k2 = ..
    k3 = ..
    k4 = ..
    
    r_n1,v_n1,theta_n1 = ..
    return r_n1,v_n1,theta_n1
    
    

def get_orbit(r0,theta0):
    bin_size = .....
    dt=.......
    
    t0=0
    t_max=.....
    
    t=[t0]
    r=[r0]
    theta=[theta0]
    
    r_n=r0
    theta_n=theta0
    v_n=0
    
    t_current=t0
    while(t_current<t_max):      
        #r_n1, v_n1, theta_n1 = Euler(dt,r_n,v_n,theta_n)  
        r_n1, v_n1, theta_n1 = rk4(dt,r_n,v_n,theta_n) 
 
        t.append(t_current)
        r.append(r_n1)
        theta.append(theta_n1)
        
        t_current += dt
        r_n = r_n1
        theta_n = theta_n1
        v_n = v_n1
    
    r1=mass1*np.array(r)/(mass1+mass2)
    r2=mass2*np.array(r)/(mass1+mass2)
    return t,r1,r2,theta

    

if __name__ == "__main__":
    mass_unit,length_unit,time_unit,G_unit = set_units(1,1,1)
    r0, theta0, period, ptheta, E_tot, mass1, mass2, mass_reduced, k_potential = get_ic(1,0,0.2)
    t,r1,r2,theta = get_orbit(r0,theta0)
    

    x1 = r1*np.sin(theta)
    y1 = r1*np.cos(theta)

    x2 = r2*np.sin(theta)
    y2 = r2*np.cos(theta)
    
    fig,ax = plt.subplots()
    ............