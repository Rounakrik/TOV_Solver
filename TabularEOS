file_path = ''
data = np.loadtxt(file_path, skiprows=7)
Solar_mass_in_km = 1.47662504

Epsilon = []
Pressure = []
NB = []

for i in range(len(data)):
    NB.append(data[i][1])
    Epsilon.append(data[i][2] * 7.4237 * (1e-19))  # Convert to Geometric Units
    Pressure.append(data[i][3] * 8.2601 * (1e-40))  # Convert to Geometric Units

Rho = np.array(NB) * 1e+54 * 1.66e-24 * 7.4237e-34  # Convert NB to Geometric Units

Pressure_function = interp1d(Epsilon, Pressure, kind='linear', fill_value='extrapolate')
Energy_Density_function = interp1d(Pressure, Epsilon, kind='linear', fill_value='extrapolate')
Rho_function = interp1d(Pressure, Rho, kind='linear', fill_value='extrapolate')

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
            (1 / 6) * (k_Epsilon_1 + 2 * k_Epsilon_2 + 2 * k_Epsilon_3 + k_Epsilon_4),
            (1 / 6) * (k_Rho_1 + 2 * k_Rho_2 + 2 * k_Rho_3 + k_Rho_4)]

# Initial conditions
Epsilon0 = 58.68e-5                               # Geometric Unit
Pressure0 = Pressure_function(Epsilon0)
Rho0 = Rho_function(Pressure0)
Radius0 = 0.01                                        # Geometric Unit
Mass0 = (4 / 3) * mth.pi * (Radius0 ** 3) * Epsilon0  # Geometric Unit
h = 0.005
Pressure_Limit = 10 ** (-14)
MassB0 = Mb0(Rho0, Mass0, Radius0) #Mass0 

MassB = [MassB0]
Mass = [Mass0]
Pressure = [Pressure0]
Epsilon = [Epsilon0]
Radius = [Radius0]
Rho = [Rho0]

# RK4 integration
while Pressure[-1] > Pressure_Limit:
    Coefficient = RK4_Slope(Mass[-1], MassB[-1], Pressure[-1], Epsilon[-1], Rho[-1], Radius[-1], h)
    Mass_i = Mass[-1] + h * Coefficient[0]
    MassB_i = MassB[-1] + h * Coefficient[1]
    Pressure_i = Pressure[-1] + h * Coefficient[2]
    Epsilon_i = Energy_Density_function(Pressure_i)
    Rho_i = Rho_function(Pressure_i)
    Radius_i = Radius[-1] + h

    MassB.append(MassB_i)
    Mass.append(Mass_i)
    Pressure.append(Pressure_i)
    Epsilon.append(Epsilon_i)
    Radius.append(Radius_i)
    Rho.append(Rho_i)

# Final results
MassB_star = MassB[-1]
Mass_star = Mass[-1]
Radius_star = Radius[-1]

print(f'Baryonic Mass of star: {MassB_star / Solar_mass_in_km} Solar masses')
print(f'Mass of star: {Mass_star / Solar_mass_in_km} Solar masses')
print(f'Radius of star: {Radius_star} km')

# Plotting results
Mass_1 = np.array(Mass)
Mass_2 = np.array(MassB)

plt.figure(figsize=(10, 6))
plt.plot(Radius, Mass_1 / Solar_mass_in_km, 'b-', label='Gravitational Mass')
plt.plot(Radius, Mass_2 / Solar_mass_in_km, 'r-', label='Baryonic Mass')
plt.xlabel('Radius (km)')
plt.ylabel('Mass (solar masses)')
plt.title('Mass vs Radius for eos_akmalpr Dataset')
# plt.savefig('Mass vs Radius for eos_akmalpr Dataset')
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(Radius, Epsilon, label='Energy Density')
plt.xlabel('Radius (km)')
plt.ylabel('Energy Density (in Geometric Unit)')
plt.title('Energy Density vs Radius for eos_akmalpr Dataset')
# plt.savefig('Energy Density vs Radius for eos_akmalpr Dataset')
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(Radius, Pressure, label='Pressure')
plt.xlabel('Radius (km)')
plt.ylabel('Pressure (in Geomtric Unit)')
plt.title('Pressure vs Radius for eos_akmalpr Dataset')
# plt.savefig('Pressure vs Radius for eos_akmalpr Dataset')
plt.show()
