from Gas import *
from Enums import *
from GasComposition import *

# for grav in range(60,171,5):
#     gg = grav / 100
#     ret = Gas_Pseudocritical(gg, 0, 0, 0, Correlation_Gas_Pseudocritical.Piper_McCain_Corredor)
#     print(f'{gg}, {ret[0]}, {ret[1]}')

# print()
# print()

# for grav in range(60,171,5):
#     gg = grav / 100
#     ret = Gas_Pseudocritical(gg, 1, 1, 1, Correlation_Gas_Pseudocritical.Piper_McCain_Corredor)
#     print(f'{gg}, {ret[0]}, {ret[1]}')

# print()
# print()

# for grav in range(60,171,5):
#     gg = grav / 100
#     ret = Gas_Pseudocritical(gg, 0, 0, 0, Correlation_Gas_Pseudocritical.Sutton_Wichert_Aziz)
#     print(f'{gg}, {ret[0]}, {ret[1]}')

# print()
# print()

# for grav in range(60,171,5):
#     gg = grav / 100
#     ret = Gas_Pseudocritical(gg, 1.0, 1.0, 1.0, Correlation_Gas_Pseudocritical.Sutton_Wichert_Aziz)
#     print(f'{gg}, {ret[0]}, {ret[1]}')


# cmp = GasComposition(145,api)

# print(ret)
        # cmp = GasComposition(119,api)
        # cmp.CO2.Frac = 0.0849
        # cmp.N2.Frac = 0.1186
        # cmp.H2S.Frac = 0.2419
        # cmp.C1.Frac = 0.3836
        # cmp.C2.Frac = 0.0629
        # cmp.C3.Frac = 0.0261
        # cmp.IC4.Frac = 0.0123
        # cmp.NC4.Frac = 0.0154
        # cmp.IC5.Frac = 0.0051
        # cmp.NC5.Frac = 0.0052
        # cmp.C6.Frac = 0.0067
        # cmp.C7Plus.Frac = 0.0373


def __GasComposition_Pseudocritical():
        api = 141.5 / 0.780 - 131.5
        cmp = GasComposition(145,api)
        cmp.CO2.Frac = 0.0642
        cmp.N2.Frac = 0.0316
        cmp.C1.Frac = 0.5026
        cmp.C2.Frac = 0.0968
        cmp.C3.Frac = 0.0620
        cmp.IC4.Frac = 0.0219
        cmp.NC4.Frac = 0.0373
        cmp.IC5.Frac = 0.0188
        cmp.NC5.Frac = 0.0185
        cmp.C6.Frac = 0.0304
        cmp.C7Plus.Frac = 0.1159
        ret = cmp.PseudoCritical()
        print(ret)

def __GasComposition_GasGravity():
    api = 141.5 / 0.758 - 131.5
    cmp = GasComposition(128,api)
    cmp.C1.Frac = 0.9712
    cmp.C2.Frac = 0.0242
    cmp.C3.Frac = 0.0031
    cmp.IC4.Frac = 0.0005
    cmp.NC4.Frac = 0.0002
    cmp.IC5.Frac = 0.0
    cmp.NC5.Frac = 0.0
    cmp.C6.Frac = 0.0002
    cmp.C7Plus.Frac = 0.0006
    ret = cmp.SpecificGravity()
    print(f'specific gravity: {ret}')

def __GasComposition_HeatingValue():
    api = 141.5 / 0.758 - 131.5
    cmp = GasComposition(103,api)
    cmp.CO2.Frac = 0.0167
    cmp.N2.Frac = 0.0032
    cmp.C1.Frac = 0.7102
    cmp.C2.Frac = 0.1574
    cmp.C3.Frac = 0.0751
    cmp.IC4.Frac = 0.0089
    cmp.NC4.Frac = 0.0194
    cmp.IC5.Frac = 0.0034
    cmp.NC5.Frac = 0.0027
    cmp.C6.Frac = 0.0027
    cmp.C7Plus.Frac = 0.0003
    ret = cmp.HeatingValue()
    print(f'Heating value: {ret}')

def __GasComposition_LiquidContent():
    api = 141.5 / 0.70 - 131.5
    cmp = GasComposition(103,api)
    cmp.CO2.Frac = 0.0167
    cmp.N2.Frac = 0.0032
    cmp.C1.Frac = 0.7102
    cmp.C2.Frac = 0.1574
    cmp.C3.Frac = 0.0751
    cmp.IC4.Frac = 0.0089
    cmp.NC4.Frac = 0.0194
    cmp.IC5.Frac = 0.0034
    cmp.NC5.Frac = 0.0027
    cmp.C6.Frac = 0.0027
    cmp.C7Plus.Frac = 0.0003
    ret = cmp.HeatingValue
    print(f'GPM: {ret}')

def __GasComposition_Viscosity():
        api = 141.5 / 0.70 - 131.5
        cmp = GasComposition(103,api)
        cmp.CO2.Frac = 0.0064
        cmp.N2.Frac = 0.0072
        cmp.C1.Frac = 0.8093
        cmp.C2.Frac = 0.099
        cmp.C3.Frac = 0.046
        cmp.IC4.Frac = 0.0076
        cmp.NC4.Frac = 0.0135
        cmp.IC5.Frac = 0.002
        cmp.NC5.Frac = 0.004
        cmp.C6.Frac = 0.0039
        cmp.C7Plus.Frac = 0.0011
        ret = cmp.Viscosity(280, 8000, 1.244)
        print(f"Viscosity: {ret}")

def __Gas_ZFact():
        for t in range(11,30):
            for p in range(2,60):
                ppr = (p*1.0) / 2.0
                tpr = (t*1.0) / 10.0
                zf = Gas_ZFactor(tpr, ppr)
                print(tpr, ppr, zf)

def __Gas_Bg():
    Tpc, Ppc = Gas_Pseudocritical(0.6698, 0.3, 1.2, 0.69, Correlation_Gas_Pseudocritical.Piper_McCain_Corredor)
    Bg = Gas_Bg(220, 6675, Tpc, Ppc, Gas_Bg_Units.rb_per_Mscf)
    print(f"Bg: {Bg}")

def __Gas_Compressibility():
    pres = 2000.0
    for t in range(11,31):
        for p in range(2, 101, 2):
            ppr = (p*1.0) / 10.0
            tpr = (t*1.0) / 10.0
            z = Gas_ZFactor(tpr, ppr)
            cg = Gas_Compressibility(tpr, pres, pres / ppr, z)
            print(tpr, ppr, cg)

def __Gas_Viscosity():
    v = Gas_Viscosity(15.67, 19.42 , 220.0)
    print(f"viscosity: {v}")

def __WetGasEquivalents():
    ggr, veq = Wet_Gas_Equivalents(0.682, 14.3027, 55.9)
    print("Total Stream:")
    print(f"Reservoir Gas Gravity: {ggr}")
    print(f"Equivalent Volume: {veq}")

    ggr, veq = Wet_Gas_Equivalents(0.679, 14.3779, 56.0, 73, 300)
    print("---------------------------------")
    print("2-Phase Separation:")
    print(f"Reservoir Gas Gravity: {ggr}")
    print(f"Equivalent Volume: {veq}")

    ggr, veq = Wet_Gas_Equivalents(0.669, 16.38941, 60.7, 73, 900, 74)
    print("---------------------------------")
    print("2-Phase Separation:")
    print(f"Reservoir Gas Gravity: {ggr}")
    print(f"Equivalent Volume: {veq}")

__WetGasEquivalents()
