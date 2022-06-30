from math import exp

Correlation_Pseudocritical = enum(Piper_McCain_Corredor = 1, Sutton_Wichert_Aziz = 2)
Correlation_ZFactor = enum(Dranchuk_Abou_Kassem = 1)
Correlation_Compressibility = enum(McCain_DAK = 1)
Correlation_Viscosity = enum(Lee_Gonzalez_Eakin = 1)
Correlation_WetGasEq = enum(Gold_McCain_Jennings = 1)
Units_Density = enum(lbm_per_cuft = 1, psi_per_ft = 2)
Units_Bg = enum(lbm_per_cuft = 1, psi_per_ft = 2)

def Gas_Pseudocritical(gg, pctH2S, pctCO2, pctN2, correlation=1):
    """Calculates pseudocrital temperature and pressure when the gas composition is not
       known.

    Arguments:
        gg:     Gas gravity (air=1)
        pctH2S: Mol % H2S (0-100)
        pctCO2: Mol % CO2 (0-100)
        pctN2:  Mol % N2 (0-100)
        correlation: correlation to use (see Enums.Correlation_Gas_Pseudocritical)
                        Piper_McCain_Corredor: 
                            Piper, L.D., McCain, W.D., Jr., and Corredor, J.H.:
                            "Compressibility Factors for Naturally Occurring Petroleum Gasses,"
                            Gas Reservoir Engineering, reprint series, 52 (1999), 186 - 200. 

                        Sutton_Wichert_Aziz:
                            Sutton, R.P.: "Compressibility Factors for High-Molecular-Weight
                            Reservoir Gasses," paper SPE 14265 presented at the 1985 SPE
                            Technical Conference and Exhibition, Las Vegas, Sept 22-25.

                            Wichert, E. and Aziz, K.: "Calculates Z's for Sour Gasses,"
                            Hyd. Proc. (May 1972) 51, 119-122.
    
    Return Value:
        Tuple(Tpc, Ppc) - Tuple containing pseudocritical temperature (deg R) and 
        pseudocritical pressure (psia)
    
    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/06/2021
    Version Notes:  Original version
    
    Known Issues: Piper, McCain, Corredor correlation - Could not find published charts to test against.
                  Tested against c# code"""

    ret = (0,0)

    if correlation == Correlation_Pseudocritical.Piper_McCain_Corredor:
        a = [1.1582E-1, -4.5820E-1, -9.0348E-1, -6.6025E-1, 7.0729E-1, -9.9397E-2]
        b = [3.8216E0, -6.5340E-2, -4.2113E-1, -9.1249E-1, 1.7438E1, -3.2191E0]
        tc = [0, 672.27, 547.54, 226.97]
        pc = [0, 1306.0, 1071.0, 493.0]
        y = [0, pctH2S/100.0, pctCO2/100.0, pctN2/100.0]

        j = a[0] + a[4] * gg + a[5] * gg ** 2
        k = b[0] + b[4] * gg + b[5] * gg ** 2

        for i in range(1,4):
            j += a[i] * y[i] * tc[i] / pc[i]
            k += b[i] * y[i] * tc[i] / pc[i] ** 0.5

        tpc = k ** 2 / j
        ppc = tpc / j

    elif correlation == Correlation_Pseudocritical.Sutton_Wichert_Aziz:
        tpc = 169.2 + 349.5 * gg - 74.0 * gg ** 2.0
        ppc = 756.8 - 131.0 * gg - 3.6 * gg ** 2.0
        a = pctH2S / 100.0 + pctCO2 / 100.0
        b = pctH2S / 100.0
        e = 120.0 * (a ** 0.9 - a ** 1.6) + 15.0 * (b ** 0.5 - b ** 4.0)
        tpc -=  e
        ppc = ppc * tpc / ((tpc + e) + b * (1 - b) * e) 
    return (tpc, ppc)


    def PseudoCritical_C7Plus(self):
        sg = 141.5 / (131.5 + self.APIC7)
        api = self.APIC7
        m = self.MWC7
        if self.TBC7 == 0.0:
            Tb = (4.5579 * m ** 0.15178 * sg ** 0.15427) ** 3
        else:
            Tb = self.TBC7
        self.TBC7 = Tb
        Tb -= 459.67
        Ppc = 2.8290406 + 0.94120109E-3 * Tb - 0.30474749E-5 * Tb ** 2.0              \
              - 0.20887611E-4 * api * Tb + 0.15184103E-8 * Tb ** 3.0                  \
              + 0.11047899E-7 * api * Tb ** 2.0 - 0.48271599E-7 * Tb * api ** 2.0     \
              + 0.13949619E-9 * (api ** 2.0) * (Tb ** 2.0)
        Tpc = 768.07121 + 1.7133693 * Tb - 0.10834003E-2 * Tb ** 2.0                  \
              - 0.89212579E-2 * api * Tb + 0.38890584E-6 * Tb ** 3.0                  \
              + 0.5309492E-5 * api * Tb ** 2.0                                        \
              + 0.327116E-7 * (api ** 2.0) * (Tb ** 2.0)
        return (Tpc, Ppc)

def Gas_ZFactor(tr, pr, correlation = 1):
    """Gas Z Factor
    
    Calculates the gas z-factor using the Dranchuk and Abou-Kassem
    correlation.  Brent's method is used to ensure convergance and
    to minimize the number of iterations.

    Arguments:
        tr: pseudoreduced temperature, dimensionless
        pr: pseudoreduced pressure, dimensionless
        correlation: Correlation_ZFactor enum
                     Dranchuk_Abou_Kassem = 1
        
    Return Value:
        gas z factor, dimensionless

    References:
        Dranchuk, P.M., and Abou-Kassem, J.H. "Calculation
        of Z Factors for Natural Gases Using Equations of
        State," Journal of Canadian Petroleum Technology 14, no 3
        (1975) 5-12.

    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/06/2021
    Version Notes:  Original version

    Known Issues: None."""

    if correlation == 1:

        def __zFunction(z, tr, pr):
            a = [0, 0.3265, -1.07, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210]
            rpr = 0.27 * pr / z / tr
            f = 1.0 + (a[1] + a[2]/tr + a[3] / tr ** 3.0 + a[4] / tr ** 4.0 + a[5] / tr ** 5.0 ) * rpr +     \
                (a[6] + a[7] / tr + a[8] / tr **2.0) * rpr ** 2.0 -                                         \
                a[9] * (a[7] / tr + a[8] / tr ** 2.0) * rpr ** 5.0 +                                        \
                a[10] * (1.0 + a[11] * rpr ** 2.0) * (rpr **2.0) / (tr ** 3.0) * exp(-a[11] * rpr ** 2.0)
            return z - f

        def __brents(f, x0, x1, tr, pr, max_iter=50, tolerance=1e-5):
            """Finds the root of a fuction f using Brent's method

            Arguments:
                f: Function to find the roots of (x where f(x)=0)
                tr: pseudoreduced temperature, dimensionless
                pr: pseudoreduced pressure, dimensionless
                x0: Lower limit for x
                x1: Upper limit for x
                max_iter: Maximum number of iterations (optional, default = 50)
                tolerance: convergence is achieved when abs(f(x)<tolerance) 
                

            Return Value:
                x: the root of the equation

            Written By:     Ray Flumerfelt (Modified from code published by Nick Ryan https://nickcdryan.com/)
            Version:        1.0
            Version Date:   03/06/2021
            Version Notes:  Original version

            Known Issues: None."""

            d=0
            fx0 = f(x0, tr, pr)
            fx1 = f(x1, tr, pr)
        
            assert (fx0 * fx1) <= 0, "Root not bracketed" 
        
            if abs(fx0) < abs(fx1):
                x0, x1 = x1, x0
                fx0, fx1 = fx1, fx0
        
            x2, fx2 = x0, fx0
        
            mflag = True
            steps_taken = 0
        
            while steps_taken < max_iter and abs(x1-x0) > tolerance:
                fx0 = f(x0, tr, pr)
                fx1 = f(x1, tr, pr)
                fx2 = f(x2, tr, pr)
        
                if fx0 != fx2 and fx1 != fx2:
                    L0 = (x0 * fx1 * fx2) / ((fx0 - fx1) * (fx0 - fx2))
                    L1 = (x1 * fx0 * fx2) / ((fx1 - fx0) * (fx1 - fx2))
                    L2 = (x2 * fx1 * fx0) / ((fx2 - fx0) * (fx2 - fx1))
                    new = L0 + L1 + L2
        
                else:
                    new = x1 - ( (fx1 * (x1 - x0)) / (fx1 - fx0) )
        
                if ((new < ((3 * x0 + x1) / 4) or new > x1) or
                    (mflag == True and (abs(new - x1)) >= (abs(x1 - x2) / 2)) or
                    (mflag == False and (abs(new - x1)) >= (abs(x2 - d) / 2)) or
                    (mflag == True and (abs(x1 - x2)) < tolerance) or
                    (mflag == False and (abs(x2 - d)) < tolerance)):
                    new = (x0 + x1) / 2
                    mflag = True
        
                else:
                    mflag = False
        
                fnew = f(new, tr, pr)
                d, x2 = x2, x1
        
                if (fx0 * fnew) < 0:
                    x1 = new
                else:
                    x0 = new
        
                if abs(fx0) < abs(fx1):
                    x0, x1 = x1, x0
        
                steps_taken += 1
        
            return x1, steps_taken

        z, steps = __brents(__zFunction, 0.25, 3.5, tr, pr)

    return z

def Gas_Bg(Tres, Pres, Tc, Pc, units=1):
    """Gas Formation Volume Factor.

    Arguments:
        Tres: Reservoir Temperature, F
        Pres: Reservoir Pressure, psia
        Tc: Critical Temperature, R
        Pc: Critical Pressure, psia
        Units: Enums.Gas_Bg_Units
               Gas_Bg_Units.rb_per_scf = 1 (reservoir bbl/scf)
               Gas_Bg_Units.rcf_per_scf = 2 (reservoir ft^3 / scf)
               Gas_Bg_Units.rb_per_Mscf = 3 (reservoir bbl/Mscf)
               Gas_Bg_Units.rcf_per_Mscf = 4 (reservoir ft^3 / Mscf)

    Return Value:
        Bg in units of Enums.Gas_Bg_Units

    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/17/2021
    Version Notes:  Original version

    Known Issues: None."""

    T = Tres + 459.67
    tr = T / Tc
    pr = Pres / Pc
    z = Gas_ZFactor(tr, pr)
    bg = z * T / Pres
    
    if units == Units_Bg.rb_per_scf:
        return bg * 0.00502
    elif units == Units_Bg.rcf_per_scf:
        return bg * 0.0282
    elif units == Units_Bg.rb_per_Mscf:
        return bg * 5.02
    else:
        return bg * 28.2

def Gas_Compressibility(Tr, Pres, Pc, z, correlation = 1):
    """Gas Coefficient of Isothermal Compressibility
       
       dz/dpr is calculated using the method presented by McCain, which
       is based on the Dranchuk Abou-Kassem equations

    Arguments:
        Tr: Pseudoreduced Temperature, dimensionless
        Pres: Reservoir Pressure, psia
        Pc: Critical Pressure, psia
        z: Z Factor, dimensionless

    Return Value:
        Gas Coefficient of Isothermal Compressibility, 1/psi 

    References:
        Dranchuk, P.M., and Abou-Kassem, J.H.,
        "Calculation of Z Factors for Natural Gasses Using Equations
        of State," Journal of Canadian Petroleum Technology 14 no. 3
        (1975), 5-12

        McCain, W.D., Jr, "The Properties of Petroleum Fluids, Third Edition" (2017), Chapter 6.

    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/17/2021
    Version Notes:  Original version

    Known Issues: None."""
    if correlation == 1:
        a = [0, 0.3265, -1.07, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210]
        pr = Pres / Pc
        rpr = 0.27 * pr / z / Tr

        dzdp = a[1] + a[2]/Tr + a[3] / Tr ** 3.0 + a[4] / Tr ** 4.0 + a[5] / Tr ** 5.0  +     \
            2.0 * rpr * (a[6] + a[7] / Tr + a[8] / Tr **2.0) -                                         \
            a[9] * 5.0 * rpr ** 4.0 * (a[7] / Tr + a[8] / Tr ** 2.0)  +                                        \
            2.0 * a[10] * rpr / Tr ** 3.0 * (1.0 + a[11] * rpr ** 2.0 - (a[11] ** 2.0) * rpr ** 4.0)  * exp(-a[11] * rpr ** 2.0)
        
        cpr = 1.0 / pr - 0.27 / z / z / Tr * dzdp / (1 + rpr / z * dzdp)
        cg = cpr / Pc
    return cg

def Gas_Density(Press, Temp, z, M, units = 1):
    """Gas Density

    Density = pM / (zRT)

    Arguments:
        Press: Pressure, psia
        Temp: Temperature, F
        z: Z Factor, dimensionless
        M: Gas Molecular Weight, lb / lb-mol (air = 29)
        units: Output units: Gas_Density_Units Enum
                lbm_per_cuft = 1
                psi_per_ft = 2

    Return Value:
        Gas Denisty in chosen units (default is lbm / ft^3) 

    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/17/2021
    Version Notes:  Original version

    Known Issues: None."""

    d = Press * M / (z * 10.732 * (Temp + 459.67))

    if units == Units_Density.psi_per_ft:
        d / 144.0
    return d

def Gas_Molecular_Weight(gg):
    """Gas Molecular Weight

    Molecular Weight = 29.0 * GasGravity

    Arguments:
        gg: Gas Gravity, air = 1

    Return Value:
        Gas Molecular Weight, lb / lb-mol

    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/18/2021
    Version Notes:  Original version

    Known Issues: None."""

    return gg * 29.0

def Gas_Viscosity(GasDensity, M, Temp, correlation = 1):
    """Gas Viscosity
    
    Calculated using the Lee, Gonzalez, and Eakin correlation

    Arguments:
        GasDensity: Gas Density, lbm / ft^3
        M: Gas molecular weight, lb / lb-mol
        Temp: Temperature, deg F
        correlation: Correlation_Viscosity enum
                     Lee_Gonzalez_Eakin = 1

    Return Value:
        Gas viscosity, cp 

    References:
        Lee, A.L., Gonzalez, M.H., and Eakin, B.E., 1966.
        "The Viscosity of Natural Gas," Journal of Petroleum Technology 18 (1966), 997-1000.

    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/17/2021
    Version Notes:  Original version

    Known Issues: None."""

    if correlation == 1:
        rho = GasDensity / 62.42796
        t = Temp + 459.67

        a = (9.379 + 0.01607 * M) * (t ** 1.5) / (209.2 + 19.26 * M + t)
        b = 3.448 + 986.4 / t + 0.01009 * M
        c = 2.447 - 0.2224 * b

        visc = a * 1E-4 * exp(b* (rho ** c))

    return visc

def Wet_Gas_Equivalents(gg, CondYield, API, SepTemp1 = 0, SepPress1 = 0,  SepTemp2 = 0, correlation = 1):
    """Calculates Reservoir Gas Gravity and Condensate Vapor Equivalent
    
    When SepTemp1 and SepPress1 are input and > 0 and SepTemp2 = 0 or not input:
        2-stage separation: gg, and CondYield is from the primary separator - does 
        not include stock tank gas
    When SepTemp1, SepPress1, and SepTemp2 are input and > 0:
        3-stage separation: gg, and CondYield is from the primary separator - does 
        not include stock tank gas
    When SepTemp1 or SepPress1 = 0 or not input:
        gg, and CondYield are assumed to be total stream gravity and yield (gas vapor is
        recovered from the stock tank and included in yield.  gg is the combined separator
        and stock tank gravity)

    Arguments:
        gg: Gas Gravity, air = 1
        CondYield: Condensate Yield, bbl/MMcf
        API: Condensate API gravity
        SepTemp1: Primary Separator Temperature, F (optional, see above)
        SepPress1: Primary Separator Pressure, psia (optional, see above)
        SepTemp2: Secondary Separator Temperature, F (optional, see above)
        correlation: 

    Return Value:
        Tuple Reservoir Gas Gravity (air=1), Condensate Vapor Equivalent (scf/stb)

    References:
        Gold, D.K., McCain, W.D., Jr., and Jennings, J.W., "An Improved Method
        for the Determination of Reservoir-Gas Specific Gravity for Retrograde
        Gases," Journal of Petroleum Technology 41 (July 1989), 747-752.

    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/21/2021
    Version Notes:  Original version

    Known Issues: None."""

    r = 1E6 / CondYield
    og = 141.5 / (131.5 + API)
    mo = 5954 / (API - 8.8)
    if SepTemp1 == 0 or SepPress1 == 0:
        veq = 133300 * og / mo
        ggr = (r * gg + 4600 * og) / (r + veq)
    elif SepTemp2 == 0:
        if correlation == Correlation_WetGasEq.Gold_McCain_Jennings:
            b0 = 635.530
            b1 = 0.361821
            b2 = 1.05435
            b3 = 5.08305
            b4 = 1.58124
            b5 = -0.791301
            veq = b0 + b1 * SepPress1 ** b2 * gg ** b3 * API ** b4 * SepTemp1 ** b5

            b1 = 1.45993
            b2 = 1.3394
            b3 = 7.09434
            b4 = 1.14356
            b5 = -0.93446
            agp = b1 * (SepPress1 - 14.65) ** b2 * gg ** b3 * API ** b4 * SepTemp1 ** b5
            ggr = (r * gg + 4600 * og + agp) / (r + veq)
    else:
        if correlation == Correlation_WetGasEq.Gold_McCain_Jennings:
            b0 = 535.916
            b1 = 2.6231
            b2 = 0.793183
            b3 = 4.6612
            b4 = 1.2094
            b5 = -0.849115
            b6 = 0.26987
            veq = b0 + b1 * SepPress1 ** b2 * gg ** b3 * API ** b4 * SepTemp1 ** b5 * SepTemp2 ** b6

            b1 = 2.99222
            b2 = 0.970497
            b3 = 6.80491
            b4 = 1.07916
            b5 = -1.19605
            b6 = 0.55367
            agp = b1 * (SepPress1 - 14.65) ** b2 * gg ** b3 * API ** b4 * SepTemp1 ** b5 * SepTemp2 ** b6
            ggr = (r * gg + 4600 * og + agp) / (r + veq)
    return (ggr, veq)


