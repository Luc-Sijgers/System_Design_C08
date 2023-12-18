# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:29:49 2023

@author: soham
"""

"""Work package 4 iteration"""
import numpy as np
import math

# Initial Constants
width = 0
diameter = 0
thickness = 0
final_mass = 0
rho = 0.284 * 0.45359237  # only for steel rn in kg/in^3
L = 0.06 * 39.37
F_ty3 = 270 * (10 ** 6) * 0.00064516129  # N/in^2 # for steel
P_ty_check = (990 * 1.5) / 4  # including safety factor and double lug
sigma_y = 0
check = 100000000000000
MS = [0]
time_t = 0
time_w = 0
time_d = 0
time_while = 0

while all(i < 0.5 for i in MS):
    check = 100000000000000
    # iterations over parameters
    for t in np.linspace(0 * 39.37, 0.02 * 39.37, 100):  # thickness
        for w in np.linspace(0 * 39.37, 0.1 * 39.37, 100):  # width
            for d in np.linspace(0.005 * 39.37, 0.05 * 39.37, 100):  # diameter
                # Check if parameters make sense and avoid divide by zero issues
                if w == d:
                    continue
                if w == 0.5 * np.sqrt(2) * d:
                    continue
                if w - 0.5 * np.sqrt(2) * d + 4 / (w - d) == 0:
                    continue
                if d > w or d == w:
                    continue
                A_av = (6 * t) / ((8 / (w - 0.5 * np.sqrt(2) * d)) + (4 / (w - d)))
                A_br = d * t
                if A_br == 0:
                    continue
                # Continue calculation
                # Correct parameters for equations MS
                x = A_av / A_br
                if x > 1.4:
                    continue
                if w / d < 1.6 or w / d > 5:
                    continue
                if (w / 2) / d < 0.5 or (w / 2) / d > 4:
                    continue
                if w < 3 * t:
                    continue
                K_ty = 0.2249 * (x ** 4) - 0.5919 * (x ** 3) + 0.1482 * (x ** 2) + 1.2507 * x - 0.0005
                P_ty = K_ty * A_br * F_ty3

                # Bending check
                I_xx = (1 / 12) * ((t / 39.37) * 1000) * (((w / 39.37) * 1000) ** 3)
                I_zz = (1 / 12) * (((t / 39.37) * 1000) ** 3) * ((w / 39.37) * 1000)
                M_z_bend = -((0.1 * 742.5) * (L / 39.37) * 1000 + 170000)
                M_x_bend = (330 * (L / 39.37) * 1000)
                sigma_y = (M_z_bend / I_zz) * (-((t / 39.37) * 1000) / 2) + (M_x_bend / I_xx) * (
                            ((w / 39.37) * 1000) / 2)

                if sigma_y > 460:
                    continue

                # Check if load is correct
                # Increase P_ty_check for MS
                if P_ty < P_ty_check + time_while:
                    continue
                mass = rho * t * (w * L + w * (d / 2) + ((math.pi * (w ** 2)) / 8) - ((math.pi * (d ** 2)) / 4))

                # Minimize mass
                if mass < check:
                    check = mass
                    final_mass = mass
                    thickness = t
                    width = w
                    diameter = d
                    # parameters for equations
                    y = A_av / A_br
                    A_t = (w - d) * t
                    e = w / 2
                    q = w / d
                    z = e / d
                    Abr = d * t
                    P_ty_test = P_ty
    print(P_ty_test)
    print('thickness (m) is: ', thickness / 39.37, 'width (m) is: ', width / 39.37, 'diameter (m) is: ',
          diameter / 39.37)
    """ print('mass (kg) is: ', final_mass)
    print('thickness (m) is: ', thickness / 39.37, 'width (m) is: ', width / 39.37, 'diameter (m) is: ', diameter / 39.37)
    print('thickness (mm) is: ', (thickness / 39.37) * 1000, 'width (cm) is: ', (width / 39.37) * 100, 'diameter (mm) is: ',
          (diameter / 39.37) * 1000) """
    # Check if margin of safety is nice
    #  equation 10 kt curves y = (A_av/A_br) [0 < y < 1.4]

    Kt1 = -0.2054 * y ** 2 + 1.4708 * y - 0.0033  # 4130 & 8630 steel 125ksi
    Kt2 = -0.631 * y ** 5 + 1.8475 * y ** 4 - 1.9182 * y ** 3 + 0.6464 * y ** 2 + 1.3237 * y + 0.0004  # 4130 & 8630 steel 150ksi
    Kt4 = -2.3872 * y ** 6 + 11.468 * y ** 5 - 19.935 * y ** 4 + 14.897 * y ** 3 - 4.9068 * y ** 2 + 1.9423 * y - 0.0004  # 4130 & 8630 steel 180ksi
    Kt5 = 0.4307 * y ** 5 - 1.6199 * y ** 4 + 2.0609 * y ** 3 - 1.3469 * y ** 2 + 1.4365 * y + 0.0004  # 356
    Kt6 = -1.9531 * y ** 6 + 8.3033 * y ** 5 - 13.371 * y ** 4 + 10.045 * y ** 3 - 3.7811 * y ** 2 + 1.6964 * y + 6 * 10 ** - 5
    Kt7 = 0.3806 * y ** 5 - 1.4623 * y ** 4 + 1.9948 * y ** 3 - 1.3698 * y ** 2 + 1.3085 * y - 8 * 10 ** -5
    Kt8 = 2.6584 * y ** 6 - 12.372 * y ** 5 + 22.117 * y ** 4 - 18.389 * y ** 3 + 6.0338 * y ** 2 + 0.6239 * y - 8E-05
    Kt9 = -2.7995 * y ** 6 + 12.87 * y ** 5 - 22.706 * y ** 4 + 19.208 * y ** 3 - 8.2669 * y ** 2 + 2.292 * y + 0.0001
    Kt11 = 0.8138 * y ** 6 - 3.9538 * y ** 5 + 6.9211 * y ** 4 - 4.7546 * y ** 3 + 0.1035 * y ** 2 + 1.3698 * y + 1E-05
    Kt12 = -1.0308 * y ** 6 + 4.7651 * y ** 5 - 8.7146 * y ** 4 + 8.3122 * y ** 3 - 4.8146 * y ** 2 + 1.9426 * y + 7E-06
    Kt13 = -0.9223 * y ** 6 + 4.9204 * y ** 5 - 10.445 * y ** 4 + 11.234 * y ** 3 - 6.493 * y ** 2 + 2.0258 * y - 2E-05
    Kt14 = 1.3021 * y ** 6 - 5.4587 * y ** 5 + 8.4936 * y ** 4 - 5.9129 * y ** 3 + 1.6035 * y ** 2 + 0.109 * y - 7E-05

    # equation 6 kt curves q = (w/d) [1 < q < 5]

    Kt1_ = -0.0015 * q ** 5 + 0.0162 * q ** 4 - 0.0621 * q ** 3 + 0.0924 * q ** 2 - 0.0686 * q + 1.0236  # 4130 & 8630 steel
    Kt2_ = -0.0015 * q ** 5 + 0.0162 * q ** 4 - 0.0621 * q ** 3 + 0.0924 * q ** 2 - 0.0686 * q + 1.0236  # 4130 & 8630 steel
    Kt3_ = -0.0015 * q ** 5 + 0.0162 * q ** 4 - 0.0621 * q ** 3 + 0.0924 * q ** 2 - 0.0686 * q + 1.0236  # 4130 & 8630 steel
    Kt4_ = -0.0005 * q ** 5 + 0.0081 * q ** 4 - 0.0522 * q ** 3 + 0.1639 * q ** 2 - 0.3324 * q + 0.8932  # AZ91C-T6 Mag, alloy sand casting, 356-T6 Aluminium Alloy Casting
    Kt5_ = 0.0003 * q ** 5 - 0.0055 * q ** 4 + 0.0403 * q ** 3 - 0.1255 * q ** 2 + 0.0674 * q + 1.0266  # 2024-T4 bar
    Kt6_ = -0.0034 * q ** 5 + 0.0491 * q ** 4 - 0.2616 * q ** 3 + 0.6192 * q ** 2 - 0.7572 * q + 1.3537  # 195T6, 220T4, and 356T6 aluminum alloy casting
    Kt7_ = 0.0017 * q ** 6 - 0.0244 * q ** 5 + 0.1232 * q ** 4 - 0.2665 * q ** 3 + 0.1938 * q ** 2 + 0.0486 * q + 0.9234  # 7075-T6 extrusion
    Kt8_ = 0.0003 * q ** 5 - 0.0055 * q ** 4 + 0.0403 * q ** 3 - 0.1255 * q ** 2 + 0.0674 * q + 1.0266  # 2024-T4 bar
    Kt9_ = 0.0017 * q ** 6 - 0.0244 * q ** 5 + 0.1232 * q ** 4 - 0.2665 * q ** 3 + 0.1938 * q ** 2 + 0.0486 * q + 0.9234  # 7075-T6 extrusion
    Kt10_ = 0.0017 * q ** 5 - 0.0285 * q ** 4 + 0.1768 * q ** 3 - 0.508 * q ** 2 + 0.5996 * q + 0.7583  # 2024-T6 Panel # 12
    Kt11_ = 0.0017 * q ** 6 - 0.0244 * q ** 5 + 0.1232 * q ** 4 - 0.2665 * q ** 3 + 0.1938 * q ** 2 + 0.0486 * q + 0.9234  # 7075-T6 extrusion
    Kt12_ = 0.0003 * q ** 5 - 0.0055 * q ** 4 + 0.0403 * q ** 3 - 0.1255 * q ** 2 + 0.0674 * q + 1.0266  # 2024-T4 bar

    # equation 7 kt curves z = (e/d) [0.5 < z < 4]

    Kbr2 = 0.0576 * z ** 3 - 0.5845 * z ** 2 + 2.4346 * z - 1.064
    Kbr3 = 0.0612 * z ** 3 - 0.621 * z ** 2 + 2.4916 * z - 1.0861
    Kbr4 = 0.0596 * z ** 3 - 0.6436 * z ** 2 + 2.5286 * z - 1.1004
    Kbr5 = 0.0746 * z ** 3 - 0.7492 * z ** 2 + 2.678 * z - 1.1562
    Kbr6 = -0.0291 * z ** 5 + 0.2929 * z ** 4 - 1.0024 * z ** 3 + 1.0035 * z ** 2 + 1.3922 * z - 0.816
    Kbr7 = -0.0192 * z ** 5 + 0.172 * z ** 4 - 0.4487 * z ** 3 - 0.1603 * z ** 2 + 2.4107 * z - 1.1137
    Kbr8 = 0.0415 * z ** 6 - 0.5046 * z ** 5 + 2.4016 * z ** 4 - 5.5347 * z ** 3 + 5.8023 * z ** 2 - 1.007 * z - 0.3632
    Kbr9 = 0.0314 * z ** 6 - 0.3693 * z ** 5 + 1.674 * z ** 4 - 3.512 * z ** 3 + 2.7785 * z ** 2 + 1.1619 * z - 0.9286
    Kbr10 = 0.0479 * z ** 5 - 0.551 * z ** 4 + 2.5122 * z ** 3 - 5.792 * z ** 2 + 7.0015 * z - 2.3904
    Kbr15 = 0.0241 * z ** 5 - 0.285 * z ** 4 + 1.3494 * z ** 3 - 3.258 * z ** 2 + 4.1168 * z - 1.3634
    Kbr20 = 0.0059 * z ** 6 - 0.0626 * z ** 5 + 0.2291 * z ** 4 - 0.2383 * z ** 3 - 0.5396 * z ** 2 + 1.6138 * z - 0.5777
    Kbr30 = 0.0031 * z ** 5 - 0.0454 * z ** 4 + 0.2715 * z ** 3 - 0.839 * z ** 2 + 1.3308 * z - 0.4397

    # Ptu calculation

    Kt = [Kt1, Kt2, Kt4, Kt5, Kt6, Kt7, Kt8, Kt9, Kt11, Kt12, Kt13, Kt14]
    F_tu = [650, 650, 650, 275, 483, 330, 573, 483, 573, 470, 573, 483]
    F_tu = [i * 10 ** 6 * 0.00064516129 for i in F_tu]  # N/in^2
    P_tu = []
    P_tu1 = []
    materials_transverse = ['4130 Steel', '4130 Steel', '4130 Steel', 'AZ91C-T6 Sand Casting', '2024-T3 Plate',
                            '220-T4 Sand casting', '7075-T6 plate', '2024-T4 Plate', '7075-T6', '2024-T42', '7075-T6',
                            '2014-T6']

    for i in range(len(Kt)):
        x = str(Kt[i] * F_tu[i] * Abr)
        P_tu.append(x + " " + materials_transverse[i])
        P_tu1.append(float(x))

    # Pbr calculation

    D_over_t_options = [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30]
    Difference = []
    Kbr = [Kbr2, Kbr3, Kbr4, Kbr5, Kbr6, Kbr7, Kbr8, Kbr9, Kbr10, Kbr15, Kbr20, Kbr30]
    P_bru = []
    Pbru1 = []

    for j in range(len(D_over_t_options)):
        diff = abs(D_over_t_options[j] - (diameter / thickness))
        Difference.append(diff)

    for i in range(len(F_tu)):
        Pbru_ = str(Kbr[Difference.index(min(Difference))] * F_tu[i] * Abr)
        Pbru1.append(float(Pbru_))
        P_bru.append(Pbru_ + " " + materials_transverse[i])

    # Pu calculation
    Pu = []
    Pu1 = []
    Kt_ = [Kt1_, Kt2_, Kt3_, Kt4_, Kt5_, Kt6_, Kt7_, Kt8_, Kt9_, Kt10_, Kt11_, Kt12_]
    for k in range(len(Kt_)):
        Pu_ = str(Kt_[k] * F_tu[k] * A_t)
        Pu.append(Pu_ + "" + materials_transverse[k])
        Pu1.append(float(Pu_))

    minPuPbr = []

    for z in range(len(Pu1)):
        if Pu1[z] < Pbru1[z]:
            minPuPbr.append(Pu1[z])
        else:
            minPuPbr.append(Pbru1[z])

    """ print('------')
    print('P_tu (N):', P_tu)  # Ptu is the ultimate transverse load for a lug only subjected to transverse loadings.
    print('------')
    print('P_bru (N):', P_bru)  # Pbru is the ultimate load due to shear-bearing loads.
    print('------')
    print('P_u (N):', Pu)  # Pu is the ultimate allowable tensile load
    print('------')
    print('minimum of Pu and Pbr (N):', minPuPbr)
    print('------') """

    # Margin of safety Work in Progress

    Fy = 329.91   # [N] Axial
    Fz = 989.73   # [N] Transverse
    Fx = Fz * 0.1  # [N] Perpendicual to both

    Rtr = [Fz / i for i in P_tu1]
    Ra = [Fy / j for j in minPuPbr]

    MS = []
    for v in range(len(Rtr)):
        MS1 = (1 / ((Ra[v] ** 1.6 + Rtr[v] ** 1.6) ** 0.625)) - 1
        MS.append(MS1)
    """print('MS is equal to:')
    print(MS)"""

    # Set new iteration boundaries
    # Make sure lower limit for thickness stays relevant
    if thickness / 39.37 == time_t:
        time_t = 0.001
    time_t += 0.001
    time_w += 0.001
    time_d += 0.001
    time_while += 25
    print(time_while)

thickness = thickness / 39.37
width = width / 39.37
diameter = diameter / 39.37

print('mass (kg) is: ', final_mass)
print('thickness (m) is: ', thickness, 'width (m) is: ', width, 'diameter (m) is: ', diameter)
print('thickness (mm) is: ', (thickness) * 1000, 'width (cm) is: ', (width) * 100, 'diameter (mm) is: ',
      (diameter) * 1000)
print(MS)

dr = 0.02  # solar panel beam thickness
# fasterner pattern
d2 = width / 5
e1 = e2 = 1.5 * d2
df = 2 * d2
length_tot = 2 * e1 + thickness + dr + thickness + 2 * e1

# bearing check
# considering coordinate axis to begin at center
# numbers go clockwise starting from bottom left
x3 = x4 = (dr / 2) + thickness + e1
x1 = x2 = -1 * x3
z2 = z3 = (width / 2) - e1
z1 = z4 = -1 * z2
fastener_coords = [(x1, z1), (x2, z2), (x3, z3), (x4, z4)]
# ____________
Fxi = Fx / 4
Fzi = Fz / 4
F_tot = np.sqrt(Fx ** 2 + Fz ** 2)
# ____________ bearing check for fastener wall
t2 = 0.002
F_ty = 460 * 1.5
sigma_br = 10000000000000000
while sigma_br > F_ty:
    sigma_br = (F_tot / (d2 * t2)) / 1000000
    if sigma_br > (F_tot / (d2 * t2)) / 1000000:
        t2 += 0.0001
# ____________ bearing check for spacecraft wall
print("t2 = ", t2)
t3 = 0.002
Fbr_sw = 345 * 1.5
while sigma_br > Fbr_sw:
    sigma_br = (F_tot / (d2 * t3)) / 1000000
    if sigma_br > (F_tot / (d2 * t2)) / 1000000:
        t3 += 0.0001

print("t3 = ", t3)


