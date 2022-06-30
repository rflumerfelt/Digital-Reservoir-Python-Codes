class GasComposition:
#region docstring
    """Gas Composition Class:
       Stores gas composition and calculates various gas properties when the gas composition
       is known.  On initialization, sets the following properties for each component, 
       based on the refence shown below (McCain):
            1.  Molecular weight (lb/lb-mol)
            2.  Pseudocritical temperature (R)
            3.  Pseudocritical pressure (psia)
        Commented out are an alternate initialization from GPA Midtream Assciation, which were
        updated in 2016.  They may be slightly better and more current than those present in
        McCain's book, but since it is likely Piper, McCain and Corredor used the values 
        published in McCain's book to develop their correlations, we have used these as a
        default.  Using different physical properties results in a difference in Tpc (+-0.3R)
        and Ppc (+-) 

    Arguments(__init__):
        MW_C7Plus:  The molecular weight of C7+ components in lb/lb-mol.  Optional, required 
                    only if C7+ mole % is > 0.
        API_C7Plus: The API gravity of C7+ components. Optional, required only if C7+ mole %
                    is > 0.
        After initialization, the mole fraction of each component present must be set.  See
        example in test cases.
    
    Reference:
        McCain, W.D., Jr, "The Properties of Petroleum Fluids, Second Edition" (1999) 492-493.
        Original Source Cited in McCain's book: "Engineering Data Book", GPSA, 1987.

        GPA Midstream Association, "Table of Physical Properties for Hydrocarbons and Other
        Compounds of Interest to the Natural Gas and Natural Gas Liquids Industries", Adopted
        as a standard 1942, Revised 2016

    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/07/2021
    Version Notes:  Original version
    
    Known Issues: The physical properties (MW, Tpc, Ppc) of individual components are set on 
                  initialization.  Various sources provide slightly different values for these
                  properties.  We chose to use the values published in McCain's book, since
                  these were almost certainly the ones used to develop the Tpc and Ppc 
                  correlations.  These properties can be overwritten after initialization.
                  See test case (Petrowiki Tpc Ppc calculations) for an example of this .

                  The calculations are very difficult to test due to the limited published
                  cases and the variability of component properties (MW, Tpc, Ppc) from different
                  sources; however, our calculations are very close with examples published by
                  McCain (but not exact).  They are reasonably close to a case published on SPE's
                  Petrowiki, but we found errors in the Petrowiki calculations.
                  More details in function documentation and test cases."""
#endregion docstring

    def __init__(self, MW_C7Plus = 0.0, API_C7Plus = 0.0):
        # Initialization using 2016 GPA Reference Table
        # for physical properties (Zsc from McCain)
        # self.C1 =  GasComponent("C1", 0.0, 16.0425, 343.01, 667.1, 0.098483, 1010.0, 0.0116, 0.0)
        # self.C2 =  GasComponent("C2", 0.0, 30.069, 549.58, 706.7, 0.077694, 1769.7, 0.0238, 0.35628)
        # self.C3 =  GasComponent("C3", 0.0, 44.0956, 665.80, 616.6, 0.072653, 2516.1, 0.0347, 0.50719)
        # self.IC4 =  GasComponent("IC4", 0.0, 58.1222, 734.06, 526.3, 0.071033, 3251.9, 0.0441, 0.56283)
        # self.NC4 =  GasComponent("NC4", 0.0, 58.1222, 765.23, 550.6, 0.071033, 3262.3, 0.0470, 0.58420)
        # self.IC5 =  GasComponent("IC5", 0.0, 72.1488, 828.63, 489.9, 0.067875, 4000.9, 0.0576, 0.62514)
        # self.NC5 =  GasComponent("NC5", 0.0, 72.1488, 845.46, 488.8, 0.069046, 4008.7, 0.0606, 0.63071)
        # self.C6 =  GasComponent("C6", 0.0, 86.1754, 913.47, 436.9, 0.0686954, 4755.9, 0.0776, 0.66406)
        # self.H2S =  GasComponent("H2S", 0.0, 34.0809, 671.58, 1305.3, 0.046125 637.1, 0.0239, 0.0)
        # self.CO2 =  GasComponent("CO2", 0.0, 44.0095, 547.43, 1070.0, 0.034257, 0.0, 0.0195, 0.0)
        # self.N2 =  GasComponent("N2", 0.0, 28.0134, 227.14, 492.5, 0.051127, 0.0, 0.00442, 0.0)
        # self.HE =  GasComponent("HE", 0.0, 4.0026, 9.35, 33.115, 0.230218, 0.0, 0.0, 0.0)

        #Initialization using reference in McCain's Book
        #for physical properties
        C7PGravity = 141.5 / (131.5 + API_C7Plus)
        self.C1 =  GasComponent("C1",0.0, 16.043, 343.00, 666.4, 0.0988, 1010, 0.0116, 0.0)
        self.C2 =  GasComponent("C2",0.0, 30.07, 549.59, 706.5, 0.0783, 1769.7, 0.0238, 0.3562)
        self.C3 =  GasComponent("C3",0.0, 44.097, 665.73, 616.0, 0.0727, 2516.2, 0.0349, 0.5070)
        self.IC4 =  GasComponent("IC4",0.0, 58.124, 734.13, 527.9, 0.0714, 3252.0, 0.0444, 0.5629)
        self.NC4 =  GasComponent("NC4",0.0, 58.124, 765.29, 550.6, 0.0703, 3262.4, 0.0471, 0.5840)
        self.IC5=  GasComponent("IC5",0.0, 72.151, 828.77, 490.4, 0.0679, 4000.9, 0.0572, 0.6247)
        self.NC5 =  GasComponent("NC5",0.0, 72.151, 845.47, 488.6, 0.0675, 4008.7, 0.0603, 0.6311)
        self.C6 =  GasComponent("C6",0.0, 86.178, 913.27, 436.9, 0.0688, 4756.0, 0.0792, 0.6638)
        self.H2S =  GasComponent("H2S",0.0, 34.076, 672.12, 1300.0, 0.0461, 637.1, 0.0239, 0.0)
        self.CO2 =  GasComponent("CO2",0.0, 44.01, 547.58, 1071.0, 0.0344, 0.0, 0.0195, 0.0)
        self.N2 =  GasComponent("N2",0.0, 28.013, 227.16, 493.1, 0.051, 0.0, 0.00442, 0.0)
        self.HE =  GasComponent("HE",0.0, 4.0026, 9.36, 33.115, 0.23, 0.0, 0.0, 0.0)
        self.C7Plus = GasComponent("C7+", 0.0, MW_C7Plus, 0.0, 0.0, 0.0, 0.0, 0.0, C7PGravity)
        self.MWC7 = MW_C7Plus
        self.APIC7 = API_C7Plus

    def GetList(self):
        """GetList: Gets the components as a list, ordered to accommodate easy calculation of Tpc and Ppc."""
        return [None, self.H2S, self.CO2, self.N2, self.C1, self.C2, self.C3, self.IC4, self.NC4, self.IC5, self.NC5, self.C6, self.C7Plus]

    def Gas_Specific_Gravity(self):
        """Calculates specific gravity of a gas mixture with known composition

        Arguments: No arguments.  The mole fractions of each component and the Molecular weight of
                    C7+ must be set before calling this function.
        
        Return Value:
            Gas Specific Gravity (Float)

        Referneces: 
                McCain, W.D., Jr, "The Properties of Petroleum Fluids, Third Edition" (2017), Chapter 5.

        Written By:     Ray Flumerfelt
        Version:        1.0
        Version Date:   03/07/2021
        Version Notes:  Original version
        
        Known Issues: None."""

        mw = self.MolecularWeight()
        return mw / 29.0

    def PseudoCritical(self):
        """Calculates pseudocrital temperature and pressure

        Arguments: No arguments.  The mole fractions of each component, Molecular weight of C7+, and API gravity of
                   C7+ must be set before calling this function.

        Return Value:
            Tuple(Tpc, Ppc) - Tuple containing pseudocritical temperature and pressure
        
        Referneces: 

               Piper, L.D., McCain, W.D., Jr., and Corredor, J.H.:
               'Compressibility Factors for Naturally Occurring Petroleum Gasses,'
               Gas Reservoir Engineering, reprint series, 52 (1999), 186 - 200. 

               Original Source:
               Piper, L.D., McCain, W.D., Jr., and Corredor, J.H.:
               'Compressibility Factors for Naturally Occurring Petroleum Gasses,' SPE 26668, Presented at the
               68th Annual Technical Conference and Exhibition of the Society of Petroleum Engineers, Houston, TX, 
               October 3-6 1993

               McCain, W.D., Jr, "The Properties of Petroleum Fluids, Third Edition" (2017), Chapter 5.
        
        Written By:     Ray Flumerfelt
        Version:        1.0
        Version Date:   03/06/2021
        Version Notes:  Original version
        
        Known Issues: None"""

        m = self.MWC7
        y = self.C7Plus.Frac
        a = [5.2073E-2, 1.0160E0, 8.6961E-1, 7.2646E-1, 8.5101E-1, 2.0818E-2, -1.506E-4]
        b = [-3.9741E-1, 1.0503E0, 9.6592E-1, 7.8569E-1, 9.8211E-1, 4.5536E-1, -3.7684E-3]
        j = a[0] + a[5] * y * m + a[6] * (y * m) ** 2.0
        k = b[0] + b[5] * y * m + b[6] * (y * m) ** 2.0
        comp = self.GetList()
        for i in range(1,4):
            j += a[i] * comp[i].Frac * comp[i].Tc / comp[i].Pc
            k += b[i] * comp[i].Frac * comp[i].Tc / comp[i].Pc ** 0.5
        for i in range(4,12):
            j += a[4] * comp[i].Frac * comp[i].Tc / comp[i].Pc
            k += b[4] * comp[i].Frac * comp[i].Tc / comp[i].Pc ** 0.5
        print(f'J = {j}, k = {k}')
        tpc = k ** 2.0 / j
        ppc = tpc / j
        return (tpc, ppc)

    def MolecularWeight(self):
        """Calculates molecular weight of a gas mixture with known composition

        Arguments: No arguments.  The mole fractions of each component and the Molecular weight of
                   C7+ must be set before calling this function.

        Return Value:
            Molecular Weight (lb/lb-mol)
        
        Referneces: 
               McCain, W.D., Jr, "The Properties of Petroleum Fluids, Third Edition" (2017), Chapter 5.
        
        Written By:     Ray Flumerfelt
        Version:        1.0
        Version Date:   03/07/2021
        Version Notes:  Original version
        
        Known Issues: None."""

        components = self.GetList()
        mw = 0
        for i in range(1,len(components)):
            mw += components[i].Frac * components[i].MW
        return mw

    def HeatingValue(self):
        """Calculates the heating value of a gas mixture with known composition

        Arguments: No arguments.  The mole fractions of each component and the Molecular weight of
                   C7+ must be set before calling this function.
        
        Return Value:
            Molecular Weight (Float)
        
        Referneces: 
               McCain, W.D., Jr, "The Properties of Petroleum Fluids, Third Edition" (2017), Chapter 5.
        
        Written By:     Ray Flumerfelt
        Version:        1.0
        Version Date:   03/07/2021
        Version Notes:  Original version
        
        Known Issues: None."""

        components = self.GetList()
        if self.C7Plus.Hv == 0:
            self.C7Plus.Hv = 53.3712 * self.C7Plus.MW + 145.418
        if self.C7Plus.Zb == 0:
            self.C7Plus.Zb = 0.00136182 * self.C7Plus.MW - 0.043665
        hv = 0
        b = 0
        for i in range(1,len(components)):
            hv += components[i].Frac * components[i].Hv
            b += components[i].Frac * components[i].Zb
        z = 1.0 - 14.65 * b ** 2
        return hv / z

    def LiquidContent(self):
        """Calculates the liquid content of a gas mixture with known composition

        Arguments: No arguments.  The mole fractions of each component, the Molecular weight 
                   and gravity of C7+ must be set before calling this function.

        Return Value:
            Liquid Content, GPM (gal/Mcf)
        
        Referneces: 
               McCain, W.D., Jr, "The Properties of Petroleum Fluids, Third Edition" (2017), Chapter 5.
        
        Written By:     Ray Flumerfelt
        Version:        1.0
        Version Date:   03/07/2021
        Version Notes:  Original version
        
        Known Issues: None."""

        components = self.GetList()
        gpm = 0
        for i in range(1,len(components)):
            if components[i].Go > 0:
                gpm += components[i].Frac * components[i].MW / components[i].Go
        return 0.3151 * gpm

    def Viscosity(self, Temp, Press, z):
        """Calculates the viscosity of a gas mixture with known composition

        Arguments: 
                    Temp: Temperature, F
                    Press: Pressure, psia
                    z: Z Factor of gas mixture, dimensionless
                    
                    The mole fractions of each component, the Molecular weight 
                    and gravity of C7+ must be set before calling this function.
        
        Return Value: Viscosity, cp 
        
        Referneces: 

               Overall Solution:
               McCain, W.D., Jr, "The Properties of Petroleum Fluids, Third Edition" (2017), Chapter 6.

               Overall Solution (original source):
               Lohrenz, J., Bray, B.G., and Clark, C.R., "Calculating Viscosity of Reservoir Fluids
               From their Composition," Journal of Petroleum Technology 16, no 10 (1964), 1171-1176.

               Calculation of Low Pressure Gas Viscosities:
               Herning, F., and Zippener, L., "Calculation of the Viscosity of Technical Gas Mixtures
               from the Viscosity of Individual Gases," Gas u. Wasserfach 79 (1936), 49-69.
               
               Calculation of Low Pressure Gas Viscosities:
               Stiel, L.I., and Thodos, G. "The Viscosity of Nonpolar Gases at Normal Pressures,"
               A.I.C.H.E.J 9 (1961), 611-615.

               Calculation of C7+ critical properties:
               Whitson, C.H., "Critical Properties Estimation from an Equation of State," paper SPE/DOE
               12634 presented at the Fourth Symposium on Enhanced Oil Recovery, Tulsa, Oklahoma,
               April 15-18, 1984.

               Convert Low Pressure Gas Viscosities to High Pressure Gas Viscosities:
               Jossi, J.A., Stiel, L.I., and Thodos, G., "The Viscosity of Pure Substances in the Dense
               Gaseous and Liquid Phases," A.I.C.H.E.J 8 (1962), 59-63.


        Return Value:
            Gas viscosity, cp
        
        Written By:     Ray Flumerfelt
        Version:        1.0
        Version Date:   03/19/2021
        Version Notes:  Original version
        
        Known Issues: None."""

        comp = self.GetList()
        t = Temp + 459.67
        api7 = self.APIC7
        rho = Press / (10.732 * z * t)
        sg7 = 141.5 / (131.5 + api7)
        kw7 = (4.5579 * self.C7Plus.MW ** 0.15178 ) * sg7 ** -0.84573
        tb7 = ((kw7*sg7) ** 3.0) - 459.67
        logpc7 = 2.8290406 + 0.94120109E-3 * tb7 - 0.30474749E-5 * tb7 ** 2    \
            - 0.20887611E-4 * api7 * tb7 + 0.15184103E-8 * tb7 ** 3.0 \
            + 0.11047899E-7 * api7 * tb7 ** 2 - 0.48271599E-7 * tb7 * (api7 ** 2) \
            + 0.13949619E-9 * (api7 ** 2) * (tb7 ** 2.0)
        self.C7Plus.Pc = 10 ** logpc7
        self.C7Plus.Tc = 768.07121 + 1.7133693 * tb7 - 0.10834003E-2 * tb7 ** 2 \
            - 0.89212579E-2 * api7 * tb7 + 0.3889058410E-6 * tb7 ** 3.0 \
            + 0.5309492E-5 * api7 * tb7 ** 2 + 0.327116E-7 * (api7 ** 2) * (tb7 ** 2)
        vm7 = 21.537 + 0.015122 * self.MWC7 - 27.656 * sg7 + 0.070615 * self.MWC7 * sg7
        self.C7Plus.Vc = vm7 / self.MWC7

        uxnum = 0
        uxden = 0
        sumv = 0
        sumtc = 0
        sumpc = 0
        summw = 0

        for i in range(1,12):
            tcj = comp[i].Tc / 1.8
            trj = t / comp[i].Tc
            pcj = comp[i].Pc / 14.69595
            prj = Press / comp[i].Pc
            Ej = tcj ** (1 / 6) * comp[i].MW ** -0.5 * pcj ** (-2 / 3)
            if trj < 1.5:
                ujx = 34E-5 / Ej * trj ** 0.94
            else:
                ujx = 17.78E-5 / Ej * (4.58 * trj - 1.67) ** 0.625
            uxnum = uxnum + comp[i].Frac * ujx * comp[i].MW ** 0.5
            uxden = uxden + comp[i].Frac * comp[i].MW ** 0.5
            sumv = sumv + comp[i].Frac * comp[i].Vc * comp[i].MW
            summw += comp[i].MW * comp[i].Frac
            sumtc += tcj * comp[i].Frac
            sumpc += pcj * comp[i].Frac
        ux = uxnum / uxden
        E = sumtc ** (1/6) * summw ** -0.5 * sumpc ** (-2/3)
        rpr = rho * sumv
        a1 = 0.1023
        a2 = .023364
        a3 = 0.058533
        a4 = -0.040758
        a5 = 0.0093324
        v = a1 + a2 * rpr + a3 * rpr ** 2 + a4 * rpr ** 3 + a5 * rpr ** 4
        v = v ** 4 - 1E-4
        v = v / E + ux

        return v

    def __str__(self):
        components = self.GetList()
        txt = "Component, Mole Fraction, Molecular Weight(lb/lb-mol), Tc (R), Pc (psia), Heating Value (BTU/scf), Summation Factor (1/psia^2), " + \
              "Liquid Gravity"
        for i in range(1,len(components)):
            txt += "\n" + components[i].__str__()
        return txt


class GasComponent:
    #region docstring
    """Gas Component Class:
       Stores physical properties for a given component

    Arguments(__init__):
        name - the component name, currently used only for printing and debugging
        frac - component mole fraction
        mw - component molecular weight (lb/lb-mol)
        tc - critical temperature (deg R)
        pc - critical pressure (psia)
        hv - gross dry heating value (btu/scf)
        zb - summation factor - used to calculate z-factor at standard conditions (1/psia^0.5)
        go - liquid specific gravity (fresh water = 1)

    Written By:     Ray Flumerfelt
    Version:        1.0
    Version Date:   03/07/2021
    Version Notes:  Original version
    
    Known Issues: None"""
#endregion docstring
    
    def __init__(self, name, frac, mw, tc, pc, vc, hv, zb, go):
        self.Name = name
        self.Frac = frac
        self.MW = mw
        self.Tc = tc
        self.Pc = pc
        self.Vc = vc
        self.Hv = hv
        self.Zb = zb
        self.Go = go

    def __str__(self):
        return f"{self.Name}, {self.Frac}, {self.MW}, {self.Tc}, {self.Pc}, {self.Vc}, {self.Hv}, {self.Zb}, {self.Go}"
