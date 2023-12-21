# Work Package 5 Iterative Code
# courtesy of Soham and Luc :D

# Imported Modules
import numpy as np


# ----------------

# GEOMETRY FUNCTIONS------------------
# MOI of cylindrical section
def I_cylindrical(R, inner):
    I_cyl = 0.25 * np.pi * (R ** 4 - inner ** 4)
    return (I_cyl)


# MOI of endcaps
def I_endcaps(R, t2):
    I_endcaps = (2 / 5) * np.pi * (R ** 4 - (R - t2) ** 4)
    return (I_endcaps)


# Mass function
def Mass_structure(R, t1, t2, rho, L):  # Mass of the whole tank
    m = np.pi * (R ** 2 - (R - t1) ** 2) * L * rho + (4 * np.pi *
                                                      (R ** 3 - (R - t2) ** 3) * rho) / 3
    return (m)


# Volume function including extra space for fins
def Vol(R, t1, t2, L):
    V = (np.pi * (R - t1) ** 2) * L + (4 * np.pi * (R - t2) ** 3) / 3 * 1.1
    return (V)


# -------------------------------------

# STRESS CHECK FUNCTIONS---------------
# Longitudinal & hoop stress----------
def Stress_long_hoop(p, R, t1, t2, sigma_yield):
    check = False
    sigma_long = (p * R) / (2 * t1)
    sigma_hoop = (p * R) / t1
    sigma_hoopt2 = (p * R) / t2
    if sigma_long < sigma_yield and sigma_hoop < sigma_yield and sigma_hoopt2 < sigma_yield:
        check = True
    return check


# Buckling Stress Check----------
def euler_buckling(rho, E, R, inner, structure_mass,
                   Fcomp):  # checking for euler stress
    check = False
    I_cyl = I_cylindrical(R, inner)
    A_cyl = np.pi * (R ** 2 - (R - t1) ** 2)
    sigma_cr = (np.pi ** 2 * E * I_cyl) / (A_cyl * L ** 2)
    sigma_xp = Fcomp / A_cyl
    if sigma_cr > sigma_xp:
        check = True
    return (check)


# Shell buckling stress----------
def Shell_buckling(P, E, R, t_1, L, nu, Fcomp):
    check = False
    A_cyl = np.pi * (R ** 2 - (R - t_1) ** 2)
    Q = (P / E) * ((R / t_1) ** 2)
    lam = np.sqrt((12 * (L ** 4) * (1 - nu ** 2)) / (np.pi ** 4 * R ** 2 * t_1 ** 2))
    k = lam + ((12 * (L ** 4) * (1 - nu ** 2)) /
               (np.pi ** 4 * R ** 2 * t_1 ** 2)) * (1 / lam)
    # Final shell buckling stress
    Sigma_cr_shell = (1.983 - 0.983 * np.exp(-23.14 * Q)) * k * (
            np.pi ** 2) * E * (1 / (12 * (1 - nu ** 2))) * ((t_1 / L) ** 2)
    sigma_xp = Fcomp / A_cyl
    if Sigma_cr_shell > sigma_xp:
        check = True
    return (check)


# Pressure Stress Checks----------
def pressure_check(sigma_yield, R, t1, p):
    check = False
    sigma_circ = p * R / t1
    sigma_lon = p * R / (2 * t1)
    if sigma_circ < sigma_yield and sigma_lon < sigma_yield:
        check = True
    return (check)


# Normal stress check----------
def normal_check(mass, L, I, A, sigma_yield):
    check = False
    sigma_tot = (6 * mass * L / 2) / I + (2 * mass) / A
    if sigma_tot < sigma_yield:
        check = True
    return (check)


# ----------------------------------

# Material Properties
p = 2310210  # N/m^2
rho = 4440  # density of material
E = 114000000000  # GPA #elastic modulus of material
nu = 0.342  # poisson ratio of material
sigma_yield = 880000000  # MPA # yield stress of material
# -----------------------------------
# Loads acting on attachments
Fz = 990  # N
Fx = 330  # N
Fy = 0  # N
F1 = 0  # N
Resultant = np.sqrt(Fz ** 2 + Fx ** 2 + Fy ** 2 + F1 ** 2)
# -----------------------------------
# Boundary Conditions and constants
L_u = 5.15 * 0.5 * 0.5  # upper limit for total tank length
R_u = 2.92 * 0.7 * 0.5  # maximum radius of tank (i suggest max diameter = half of structure diameter)
min_vol = 0.4858  # m^3 #minimum volume required for the fuel
mass_fuel = 699.9 + 424  # KG #mass of the fuel
# -----------------------------------
Mass_attachment = 0.1476  # KG
attachments = 4  # number of attachment

# Iterations
min_check = 10000000000000000000000

# ITERATION CODE
for t1 in np.arange(0.001, 0.008, 0.0002):
    for t2 in np.arange(0.0005, 0.008, 0.0002):
        for R in np.arange(0.1, R_u, 0.01):
            for L in np.arange(0.001, L_u, 0.005):
                # Checks for nonsensical values
                L_tot = L + 2 * R
                if L_tot > L_u:  # length check
                    continue
                mass = Mass_structure(R, t1, t2, rho, L)
                mass_tot = mass + mass_fuel
                I_tot = I_cylindrical(R, R - t1) + I_endcaps(R, t2)
                A = np.pi * (R ** 2 - (R - t1) ** 2)
                F_comp = mass_tot * 6.4 * 9.81
                vol = Vol(R, t1, t2, L)
                if vol < min_vol:  # volume check
                    continue

                # Stress checks
                check_1 = euler_buckling(rho, E, R, R - t1, mass, F_comp)
                check_2 = Shell_buckling(p, E, R, t1, L, nu, F_comp)
                check_4 = normal_check(mass_tot, L, I_tot, A, sigma_yield)
                check_5 = Stress_long_hoop(p, R, t1, t2, sigma_yield)
                if check_1 == True and check_2 == True and check_4 == True and check_5 == True:
                    if vol * mass < min_check:  # min mass and vol check
                        min_check = vol * mass
                        Volume, t, t_end, Len, Radius, Mass = vol, t1, t2, L, R, mass

print(t, t_end, Radius, Len, Volume, Mass)
