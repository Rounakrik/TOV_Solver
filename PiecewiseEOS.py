file_path = 
Gamma_array = np.loadtxt(file_path, skiprows = 2, usecols = 0, max_rows = 7)
K0 = 6.8011 *10 **(-9)
LogRho = np.loadtxt(file_path, skiprows = 12, usecols = 0, max_rows = 6)
Rho_array = np.array([10**i for i in LogRho])  #in g cm-3


K = [K0] # K for cgs outputs
for i in range(1,len(Rho_array)+1):
    Ki = K[i-1]*((Rho_array[i-1])**(Gamma_array[i-1]-Gamma_array[i]))
    K.append(Ki)
    
for i in range(0,len(K)):
    K[i] *= (7.4237e-19)**(1-Gamma_array[i]) #for geometric outputs
    
Rho_array = Rho_array*(7.4237e-19) #in km-2

def dmdr(Epsilon, Radius):
    return 4 * mth.pi * Epsilon * (Radius ** 2)

def dpdr(Epsilon, Pressure, Mass, Radius):
    return (-(Epsilon + Pressure) * (Mass + (4 * mth.pi * (Radius ** 3) * Pressure))) / (Radius * (Radius - (2 * Mass)))

def dMdr(Rho, Mass, Radius):
    return 4 * mth.pi * Rho * (Radius ** 2) * ((1 - 2 * Mass / Radius) ** (-1 / 2))

def Mb0(rho0, m0, r0):
    M = [m0]
    r = r0
    i = 0
    h = 0.005  # Step size
    while r <= r0:
        K_Mass_B_1 = dMdr(rho0, m0, r)
        K_Mass_B_2 = dMdr(rho0, m0, r + 0.5 * h)
        K_Mass_B_3 = dMdr(rho0, m0, r + 0.5 * h)
        K_Mass_B_4 = dMdr(rho0, m0, r + h)
        m = M[i] + ((1 / 6) * (K_Mass_B_1 + 2 * K_Mass_B_2 + 2 * K_Mass_B_3 + K_Mass_B_4))
        M.append(m)
        r += h
        i += 1
    return M[-1]

def Pressure_function(Rho): #in km-2
    for i in range(len(Rho_array) - 1): 
        if Rho>0 and Rho<=Rho_array[0]:
            return K[0] * (Rho ** Gamma_array[0])
        elif Rho_array[i] < Rho <= Rho_array[i + 1]:
            return K[i+1] * (Rho ** Gamma_array[i+1])
        elif Rho > Rho_array[-1]:
            return K[-1] * (Rho ** Gamma_array[-1])

Pressure_array = np.array([Pressure_function(i) for i in Rho_array])

def Rho_function(Pressure): #in km-2
    for i in range(len(Pressure_array)-1):
        if Pressure <=Pressure_array[0]:
            return (Pressure/K[0])**(1/Gamma_array[0])        
        elif Pressure_array[i] < Pressure <=Pressure_array[i+1]:
            return (Pressure/K[i+1])**(1/Gamma_array[i+1])
        elif Pressure > Pressure_array[-1]:
            return (Pressure/K[-1])**(1/Gamma_array[-1])
        
A_array = [0] #ai

for i in range(1,len(K)):
    ai = A_array[i-1] + ((Pressure_array[i-1]/Rho_array[i-1])*((1/(Gamma_array[i-1]-1))-(1/(Gamma_array[i]-1))))
    A_array.append(ai)

def Energy_Density_function(Pressure):
    for i in range(len(Pressure_array)-1):
        if Pressure <= Pressure_array[0]:
            return ((1+A_array[0]) * Rho_function(Pressure)+((K[0]*(Rho_function(Pressure)**(Gamma_array[0]))))/(Gamma_array[0]-1))
        elif Pressure_array[i] < Pressure <= Pressure_array[i+1]:
            return ((1+A_array[i+1])*Rho_function(Pressure))+((K[i+1]*(Rho_function(Pressure)**(Gamma_array[i+1])))/(Gamma_array[i+1]-1))
        elif Pressure >= Pressure_array[-1]:
            return ((1+A_array[-1])*Rho_function(Pressure))+((K[-1]*(Rho_function(Pressure)**(Gamma_array[-1])))/(Gamma_array[-1]-1))
    

def RK4_Slope(Mass, MassB, Pressure, Epsilon, Rho, Radius, h):
    k_Mass_1 = dmdr(Epsilon, Radius)
    k_Mass_B_1 = dMdr(Rho, Mass, Radius)
    k_Pressure_1 = dpdr(Epsilon, Pressure, Mass, Radius)
    k_Epsilon_1 = Energy_Density_function(Pressure)
    k_Rho_1 = Rho_function(Pressure)

    k_Mass_2 = dmdr(Epsilon + 0.5 * h * k_Epsilon_1, Radius + 0.5 * h)
    k_Mass_B_2 = dMdr(Rho + 0.5 * h * k_Rho_1, Mass + 0.5 * h * k_Mass_1, Radius + 0.5 * h)
    k_Pressure_2 = dpdr(Epsilon + 0.5 * h * k_Epsilon_1, Pressure + 0.5 * h * k_Pressure_1, Mass + 0.5 * h * k_Mass_1, Radius + 0.5 * h)
    k_Epsilon_2 = Energy_Density_function(Pressure + 0.5 * h * k_Pressure_1)
    k_Rho_2 = Rho_function(Pressure + 0.5 * h * k_Pressure_1)

    k_Mass_3 = dmdr(Epsilon + 0.5 * h * k_Epsilon_2, Radius + 0.5 * h)
    k_Mass_B_3 = dMdr(Rho + 0.5 * h * k_Rho_2, Mass + 0.5 * h * k_Mass_2, Radius + 0.5 * h)
    k_Pressure_3 = dpdr(Epsilon + 0.5 * h * k_Epsilon_2, Pressure + 0.5 * h * k_Pressure_2, Mass + 0.5 * h * k_Mass_2, Radius + 0.5 * h)
    k_Epsilon_3 = Energy_Density_function(Pressure + 0.5 * h * k_Pressure_2)
    k_Rho_3 = Rho_function(Pressure + 0.5 * h * k_Pressure_2)

    k_Mass_4 = dmdr(Epsilon + h * k_Epsilon_3, Radius + h)
    k_Mass_B_4 = dMdr(Rho + h * k_Rho_3, Mass + h * k_Mass_3, Radius + h)
    k_Pressure_4 = dpdr(Epsilon + h * k_Epsilon_3, Pressure + h * k_Pressure_3, Mass + h * k_Mass_3, Radius + h)
    k_Epsilon_4 = Energy_Density_function(Pressure + h * k_Pressure_3)
    k_Rho_4 = Rho_function(Pressure + h * k_Pressure_3)

    return [(1 / 6) * (k_Mass_1 + 2 * k_Mass_2 + 2 * k_Mass_3 + k_Mass_4),
            (1 / 6) * (k_Mass_B_1 + 2 * k_Mass_B_2 + 2 * k_Mass_B_3 + k_Mass_B_4),
            (1 / 6) * (k_Pressure_1 + 2 * k_Pressure_2 + 2 * k_Pressure_3 + k_Pressure_4),
           ]

print(f'Pressure points: {Pressure_array}')
print(f'Rest Mass Density in Geometric Unit: {Rho_array}')
print(f'Values of K: {K}')
print(f'Values of Gamma: {Gamma_array}')
print(f'Values of A: {A_array}')
