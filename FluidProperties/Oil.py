from math import *

def enum(**named_values):
    return type('Enum', (), named_values)

Correlation_Bubble_Point = enum(Valko_McCain = 1, Standing = 2, AlMarhoun = 3)
Correlation_Solution_GOR = enum(Velarde_McCain_Blasingame = 1, Standing = 2, AlMarhoun = 3)
Correlation_Oil_Compressibility = enum(Spivey_Valko_McCain = 1)
Correlation_Oil_Density = enum(McCain_Hill = 1)
Correlation_Oil_FVF = enum(McCain = 1)
Correlation_Oil_Viscosity = enum(Beggs_Robinson_Petrosky = 1)
Compressibility_Application = enum(pta_modeling = 1, pvt = 2, matbal = 3)

def Bubble_Point(Rsb, API, gg, Temp, correlation = 1, StandingCoeffs = [18.2, 1.4, 0.83, 0.00091, 0.0125]):
    """Calculates oil bubble point pressure.

    Arguments:
        Rsb: Bubble point solution gas/oil ratio, scf/stb
             Valko & Al-Martoun - use separator gas/oil ratio
             Standing - use total gas/oil ratio
        API: Oil API gravity
        gg:  Solution gas gravity, air = 1
        Temp: Temperature, F
        correlation: correlation to use (Enum Correlation_Bubble_Point)
                     1 = Valko and McCain
                     2 = Standing
                     3 = Al-Marhoun
    
    Return Value:
        Bubble point pressure, psia

    References:
        Valko, P.P., and McCain, W. D., Jr, "Reservoir Oil Bubblepoint Pressures
        Revisited; Solution Gas/Oil Ratios and Surface Gas Specific Gravities,"
        Journal of Petroleum Science and Engineering, 37(3): 153-169 (2003).

        Standing, M.B., "Volumetric and Phase Behavior of Oil Field Hydrocarbon
        Systems, SPE Dallas (1977).

        Al-Marhoun, M.A., "PVT Correlations for Middle East Crude Oils", Journal
        of Petroleum Technology (May 1988) 40, 650-666.
    
    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/21/2021
    Version Notes:  Original version
    
    Known Issues: None"""

    if correlation == Correlation_Bubble_Point.Valko_McCain:  
        v = log(Rsb)
        z1 = -5.48 - 0.0378 * v + 0.281 * v ** 2 - 0.0206 * v ** 3
        v = API
        z2 = 1.27 - 0.0449 * v + 4.36E-4 * v ** 2 - 4.76E-6 * v ** 3
        v = gg
        z3 = 4.51 - 10.84 * v + 8.39 * v ** 2 - 2.34 * v ** 3
        v = Temp
        z4 = -0.7835 + 6.23E-3 * v - 1.22E-5 * v ** 2 + 1.03E-8 * v ** 3
        z = z1 + z2 + z3 + z4
        pb = exp(7.475 + 0.713 * z + 0.0075 * z ** 2)
    
    elif correlation == Correlation_Bubble_Point.Standing:
        a = StandingCoeffs
        cn = (Rsb / gg) ** a[2] * 10 ** (a[3] * Temp - a[4] * API)
        pb = a[0] * (cn - a[1])

    elif correlation == Correlation_Bubble_Point.AlMarhoun:
        og = 141.5 / (131.5 + API)
        pb = 5.38088E-3 * Rsb ** 0.715082 * gg ** -1.877840 * og ** 3.1437 * (Temp + 459.67) ** 1.32657
    else:
        pb = None
    
    return pb

def Solution_GOR(Pb, Press, API, gg, Temp, Rsb = 0, correlation = 1, StandingCoeffs = [18.2, 1.4, 0.83, 0.00091, 0.0125]):
    p = Press
    if p > Pb:
        p = Pb
    if correlation == Correlation_Solution_GOR.Velarde_McCain_Blasingame:
        a1 = 9.73E-7 * gg ** 1.672608 * API ** 0.92987 * Temp ** 0.247235 * (Pb-14.7) ** 1.056052
        a2 = 0.022339 * gg ** -1.00475 * API ** 0.337711 * Temp ** 0.132795 * (Pb-14.7) ** 0.302065
        a3 = 0.725167 * gg ** -1.48548 * API ** -0.164741 * Temp ** -0.09133 * (Pb-14.7) ** 0.047094
        pr = (p - 14.7) / (Pb - 14.7)
        rsr = a1 * pr ** a2 + (1 - a1) * pr ** a3
        rs = rsr * Rsb
    elif correlation == Correlation_Solution_GOR.Standing:
        a = StandingCoeffs
        rs = ((p / a[0] + a[1]) * 10 ** (a[4] * API) / 10 ** (a[3] * Temp)) ** (1/a[2]) * gg
    elif correlation == Correlation_Solution_GOR.AlMarhoun:
        og = 141.5 / (131.5 + API)
        a1 = 5.38088E-3
        a2 = 0.715082
        a3 = -1.87784
        a4 = 3.1437
        a5 = 1.32657
        rs = (p / a1 / gg ** a3 / og ** a4 / (Temp + 459.67) ** a5) ** (1/a2)
    return rs

def Oil_Viscosity(Pb, Press, Rs, API, Temp, correlation = Correlation_Oil_Viscosity.Beggs_Robinson_Petrosky):
    a = 10.715 * (Rs + 100) ** -0.515
    b = 5.44 * (Rs + 150) ** -0.338
    c = 10 ** (3.0324 - 0.02023 * API) * Temp ** -1.163
    mod = 10 ** c - 1
    mo = a * mod ** b
    if Press > Pb:
        a = -1.0146 + 1.3322 * log10(mo) - 0.4876 * (log10(mo)) ** 2 - 1.15036 * (log10(mo)) ** 3
        mo = mo + 0.0013449 * (Press - Pb) * 10 ** a
    return mo

def Oil_Reservoir_Gas_Gravity(Press, Rsb, API, ggsep1, Temp):
    a = -208.0797 / Press + 22885 / Press ** 2 - 0.000063641 * Press + 3.38346 / Temp ** 0.5 - 0.000992 * Temp \
        - 0.000081147 * Rsb - 0.001956 * API + 1.081956 / ggsep1 + 0.394035 * ggsep1 ** 2
    return 1 / a

def Oil_Compressibility_Undersaturated(Pb, Press, Rsb, API, gg, Temp, correlation = 1, application = 1, Pinit = 0):
    if correlation == Correlation_Oil_Compressibility.Spivey_Valko_McCain:
        def cofb(Pb, Press, Rsb, API, gg, Temp):
            v = log(API)
            z1 = 3.011 - 2.6254 * v + 0.497 * v ** 2
            v = log(gg)
            z2 = -0.0835 - 0.259 * v + 0.382 * v ** 2
            v = log(Pb)
            z3 = 3.51 - 0.0289 * v - 0.0584 * v ** 2
            v = log(Press/Pb)
            z4 = 0.327 - 0.608 * v + 0.0911 * v ** 2
            v = log(Rsb)
            z5 = -1.918 - 0.642 * v + 0.154 * v ** 2
            v = log(Temp)
            z6 = 2.52 - 2.73 * v + 0.429 * v ** 2
            z = z1 + z2 + z3 + z4 + z5 + z6
            return exp(2.434 + 0.475 * z + 0.048 * z ** 2 - log(1E6)), z
        if Press >= Pb:
            if application == Compressibility_Application.pta_modeling:
                cp, z = cofb(Pb, Press, Rsb, API, gg, Temp)
                dzdp = (-0.608 + 0.1822 * log(Press / Pb)) / Press
                dcdp = cp * (0.475 + 0.096 * z) * dzdp
                return cp + (Press - Pb) * dcdp
            elif application == Compressibility_Application.pvt:
                return cofb(Pb, Press, Rsb, API, gg, Temp)[0]
            elif application == Compressibility_Application.matbal:
                cp = cofb(Pb, Press, Rsb, API, gg, Temp)[0]
                cpi = cofb(Pb, Pinit, Rsb, API, gg, Temp)[0]
                return ((Press - Pb) * cp - (Pinit - Pb) * cpi) / (Press - Pinit)

def Oil_Compressibility_Saturated(Pb, Press, Rsb, API, gg, ggsep1, Bo, Bg, Temp):
    og = 141.5 / (131.5 + API)
    a1 = 9.73E-7 * ggsep1 ** 1.672608 * API ** 0.92987 * Temp ** 0.247235 * (Pb-14.7) ** 1.056052
    a2 = 0.022339 * ggsep1 ** -1.00475 * API ** 0.337711 * Temp ** 0.132795 * (Pb-14.7) ** 0.302065
    a3 = 0.725167 * ggsep1 ** -1.48548 * API ** -0.164741 * Temp ** -0.09133 * (Pb-14.7) ** 0.047094
    pr = (Press - 14.7) / (Pb - 14.7)
    rsr = a1 * pr ** a2 + (1 - a1) * pr ** a3
    Rs = rsr * Rsb
    drsdp = Rsb * (a1 * a2 * pr ** (a2 - 1) + (1 - a1) * a3 * pr ** (a3 - 1)) / (Pb - 14.7)

    rpo = 52.8 - 0.01 * Rsb
    ra = 0
    rpold = 10000
    count = 0
    while abs(rpo-rpold) > 0.001 or count > 15:
        count = count + 1
        rpold = rpo
        ra = -49.893 + 85.0149 * ggsep1 - 3.70373 * ggsep1 * rpo + 0.0479818 * ggsep1 * rpo ** 2 + 2.98914 * rpo - 0.0356888 * rpo ** 2
        rpo = (Rsb * gg + 4600 * og) / (73.71 + Rsb * gg / ra)
    drp = (0.167 + 16.181 * 10 ** (-0.0425 * rpo)) * (Press / 1000) - 0.01 * (0.299 + 263 * (10 ** (-0.0603 * rpo))) * (Press / 1000) ** 2
    rbs = rpo + drp
    drt = (0.00302 + 1.505 * rbs ** -0.951) * (Temp - 60) ** 0.938 - (0.0216 - 0.0233 * 10 ** (-0.0161 * rbs)) * (Temp - 60) ** 0.475
    ror = rbs - drt
    rbs = rpo + drp
    drodp = gg * drsdp * (73.71 - 4600 * og / ra) / (73.71 + Rs * gg / ra) ** 2
    drpdp = 1E-3 * (0.167 + 16.181 * 10 ** (-0.0425 * rpo)) - \
            1E-3 * (1.5835 * 10 ** (-0.0425 * rpo) * Press * drodp) - \
            1E-8 * (0.598 * Press + 526 * Press * 10 ** (-0.0603 * rpo)) + \
            1E-8 * (36.52 * Press ** 2 * 10 ** (-0.0603 * rpo) * drodp)
    drtdp = (-1.4313 * rbs ** -1.9551 * (Temp - 60) ** 0.938 - 0.0008638 * 10 ** (-0.0161 * rbs) * (Temp - 60) ** 0.475) * (drodp + drpdp)
    drordp = drodp + drpdp + drtdp
    dbodp = 1 / ror ** 2 * (0.01357 * gg * ror * drsdp - (og + 0.01357 * Rs * gg) * drordp)
    co = -1 / Bo * (dbodp - Bg * drsdp)
    return co

def Oil_Density(Pb, Press, Rsb, API, ggsep1, ggtot, Temp, correlation = Correlation_Oil_Density.McCain_Hill):
    p = Press
    if p > Pb:
        p = Pb
    og = 141.5 / (131.5 + API)
    rpo = 52.8 - 0.01 * Rsb
    ra = 0
    rpold = 10000
    count = 0
    while abs(rpo-rpold) > 0.001 or count > 15:
        count = count + 1
        rpold = rpo
        ra = -49.893 + 85.0149 * ggsep1 - 3.70373 * ggsep1 * rpo + 0.0479818 * ggsep1 * rpo ** 2 + 2.98914 * rpo - 0.0356888 * rpo ** 2
        rpo = (Rsb * ggtot + 4600 * og) / (73.71 + Rsb * ggtot / ra)
    drp = (0.167 + 16.181 * 10 ** (-0.0425 * rpo)) * (p / 1000) - 0.01 * (0.299 + 263 * (10 ** (-0.0603 * rpo))) * (p / 1000) ** 2
    rbs = rpo + drp
    drt = (0.00302 + 1.505 * rbs ** -0.951) * (Temp - 60) ** 0.938 - (0.0216 - 0.0233 * 10 ** (-0.0161 * rbs)) * (Temp - 60) ** 0.475
    rob = rbs - drt
    ro = rob
    if Press > Pb:
        co = Oil_Compressibility(Pb, Press, Rsb, API, ggsep1, Temp, correlation = Correlation_Oil_Compressibility.Spivey_Valko_McCain, \
             application=Compressibility_Application.pvt)
        ro = rob * exp(co * (Press - Pb))
    return ro

def Oil_FVF(Pb, Press, Rs, API, gg, Temp = 0, correlation = Correlation_Oil_FVF.McCain, Co = 0, Density = 0):
    rsto = 141.5 / (131.5 + API) * 62.37
    bo = (rsto + 0.01357 * Rs * gg)  / Density
    return bo





