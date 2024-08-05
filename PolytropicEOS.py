import numpy as np
import math as mth
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['text.usetex'] = True
rcParams["font.family"] = "Times New Roman"

# G = 1
# c = 1
Solar_mass_in_km = 1.47662504

def Energy_Density_function(Pressure):
    K = 100
    Gamma = 2
    return (Pressure/K) ** (1/Gamma)

def Pressure_function(Epsilon):
    K = 100
    Gamma = 2
    return K*((Epsilon)**(Gamma))

def dmdr(Epsilon, Radius):
    return 4*mth.pi*Epsilon*(Radius**2)

def dpdr(Epsilon, Pressure, Mass, Radius):
    return (-(Epsilon + Pressure )*(Mass + ((4*mth.pi*(Radius**3))*Pressure)))/(Radius*(Radius-(2*Mass)))


Epsilon0 = 6*58.68e-5                        # Geometric Unit
Radius0 = 0.1                                # Geometric Unit
Mass0 = (4/3)*(mth.pi)*(Radius0**3)*Epsilon0 # Geometric Unit
Pressure0 = Pressure_function(Epsilon0)
h = 0.005
Pressure_Limit = 10**(-14)


def RK4_Slope(Mass, Epsilon, Pressure, Radius, h):
    
    k_Mass_1     = dmdr(Epsilon, Radius)
    k_Pressure_1 = dpdr(Epsilon, Pressure, Mass, Radius)
    k_Epsilon_1  = Energy_Density_function(Pressure)
   
    k_Mass_2     = dmdr(Epsilon + 0.5 * h * k_Epsilon_1, Radius + 0.5 * h)
    k_Pressure_2 = dpdr(Epsilon + 0.5 * h * k_Epsilon_1, Pressure + 0.5 * h * k_Pressure_1, Mass + 0.5 * h * k_Mass_1, Radius + 0.5 * h)
    k_Epsilon_2  = Energy_Density_function(Pressure + 0.5 * h * k_Pressure_1)
    
    k_Mass_3     = dmdr(Epsilon + 0.5 * h * k_Epsilon_2, Radius + 0.5 * h)
    k_Pressure_3 = dpdr(Epsilon + 0.5 * h * k_Epsilon_2, Pressure + 0.5 * h * k_Pressure_2, Mass + 0.5 * h * k_Mass_2, Radius + 0.5 * h)
    k_Epsilon_3  = Energy_Density_function(Pressure + 0.5 * h * k_Pressure_2)
    
    k_Mass_4     = dmdr(Epsilon + h * k_Epsilon_3, Radius + h)
    k_Pressure_4 = dpdr(Epsilon + h * k_Epsilon_3, Pressure + h * k_Pressure_3, Mass + h * k_Mass_3, Radius + h)
    k_Epsilon_4  = Energy_Density_function(Pressure + h * k_Pressure_3)
    
    return [(1/6)*(k_Mass_1+(2*k_Mass_2)+(2*k_Mass_3)+k_Mass_4), (1/6)*(k_Epsilon_1+(2*k_Epsilon_2)+(2*k_Epsilon_3)+k_Epsilon_4), (1/6)*(k_Pressure_1+(2*k_Pressure_2)+(2*k_Pressure_3)+k_Pressure_4)]

Mass = [Mass0]
Pressure = [Pressure0]
Epsilon = [Epsilon0]
Radius = [Radius0]


while Pressure[-1] > Pressure_Limit:
    Coefficient = RK4_Slope(Mass[-1], Epsilon[-1], Pressure[-1], Radius[-1], h)
    Mass_i = Mass[-1] + h * Coefficient[0]
    Pressure_i = Pressure[-1] + h * Coefficient[2]
    Epsilon_i = Energy_Density_function(Pressure_i)
    Radius_i = Radius[-1] + h
    Mass.append(Mass_i)
    Pressure.append(Pressure_i)
    Epsilon.append(Epsilon_i)
    Radius.append(Radius_i)

print(f'Mass of star: {Mass[-1] / Solar_mass_in_km}')
print(f'Radius of star: {Radius[-1]}')


Mass = np.array(Mass)

plt.figure(figsize=(10, 5))
plt.plot(Radius, Mass / 1.476)
plt.xlabel(r'$Radius \rho (km)$')
#plt.ylabel('Mass (solar masses)')
#plt.legend()
plt.show()


plt.figure(figsize=(10, 5))
plt.plot(Radius, Pressure, label='Mass')
plt.xlabel('Radius (km)')
plt.ylabel('Pressure')
plt.title('Pressure vs Radius ')
plt.show()


plt.figure(figsize=(10, 5))
plt.plot(Radius, Epsilon, label='Energy Density')
plt.xlabel('Radius (km)')
plt.ylabel('Mass (solar masses)')
plt.title('Energy Density vs Radius ')
plt.show()



