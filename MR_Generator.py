Solar_mass_in_km = 1.47662504

final_masses_s = []
final_radii_s = []
final_massb_s = []


Epsilon0_s = np.linspace(0.4, 4, 100) * 58.68e-5
Pressure0_s = np.array([Pressure_function(i) for i in Epsilon0_s])
Rho0_s = np.array([Rho_function(i) for i in Pressure0_s])

for n in Epsilon0_s:
    
    
    Epsilon0 = n                                          # Geometric Unit
    Pressure0 = Pressure_function(Epsilon0)
    Rho0 = Rho_function(Pressure0)
    Radius0 =                                             # Put the desired central value of radius in Geometric Unit
    Mass0 = (4 / 3) * mth.pi * (Radius0 ** 3) * Epsilon0  # Geometric Unit
    h = 0.005                                             # Step size
    Pressure_Limit = 10 ** (-14)
    MassB0 = Mb0(Rho0, Mass0, Radius0)


    MassB = [MassB0]
    Mass = [Mass0]
    Pressure = [Pressure0]
    Epsilon = [Epsilon0]
    Radius = [Radius0]
    Rho = [Rho0]

    
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

    final_massb_s.append(MassB[-1] / Solar_mass_in_km)
    final_masses_s.append(Mass[-1] / Solar_mass_in_km)
    final_radii_s.append(Radius[-1])


plt.figure(figsize=(10, 6))
plt.plot(final_radii_s, final_massb_s, 'r-', label='Baryonic Mass')
plt.plot(final_radii_s, final_masses_s, 'b-', label='Gravitational Mass')
plt.xlabel('R [km]', fontsize=15)
plt.ylabel('M $[M_{\odot}]$', fontsize=15)
plt.title('M-R Curve')
plt.legend()
plt.show()
