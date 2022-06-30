from Oil import *

def Bubble_Point_Tests():
    #From McCain 3rd Edition
    pb = Bubble_Point(981, 41.5, 0.782, 184, Correlation_Bubble_Point.Valko_McCain)
    assert(abs(pb - 2993) < 1.0)
    print (f"Pb: {pb}")

    #From McCain 2nd Edition
    pb = Bubble_Point(768, 40.7, 0.786, 220, Correlation_Bubble_Point.Standing)
    assert(abs(pb - 2685) < 1.0)
    print (f"Pb: {pb}")

    #No data to validate Al_Marhoun bubble point
    pb = Bubble_Point(768, 40.7, 0.786, 220, Correlation_Bubble_Point.AlMarhoun)
    print (f"Pb: {pb}")

def Solution_GOR_Tests():
    #From McCain 3rd Edition
    rs = Solution_GOR(2993, 1200, 41.5, 0.782, 184, 981, Correlation_Solution_GOR.Velarde_McCain_Blasingame)
    assert(abs(rs - 441) < 1.0)
    print (f"Rs: {rs}")

    rs = Solution_GOR(1591.7, 614.7, 25.5, 0.76464, 272, 185 )
    #assert(abs(rs - 441) < 1.0)
    print (f"Rs example: {rs}")
    #From McCain 2nd Edition
    rs = Solution_GOR(2685.77, 2685, 40.7, 0.786, 220, 981, Correlation_Solution_GOR.Standing)
    assert(abs(rs - 768) < 1.0)
    print (f"Rs: {rs}")

    #No data to validate Al_Marhoun bubble point
    rs = Solution_GOR(3017.58, 3017.58, 40.7, 0.786, 220, 981, Correlation_Solution_GOR.AlMarhoun)
    assert(abs(rs - 768) < 1.0)
    print (f"Rs: {rs}")

def Oil_Compressibility_Undersaturated_Tests():
    co = Oil_Compressibility_Undersaturated(2993, 4000, 981, 41.5, .782, 184, \
         correlation = Correlation_Oil_Compressibility.Spivey_Valko_McCain, \
         application = Compressibility_Application.pvt)
    assert(abs(co * 1E8 - 1668) < 1.0)
    print (f"Co, Spivey - pvt: {co}")

    co = Oil_Compressibility_Undersaturated(2993, 4000, 981, 41.5, .782, 184, \
         correlation = Correlation_Oil_Compressibility.Spivey_Valko_McCain, \
         application = Compressibility_Application.matbal, Pinit = 5000)
    assert(abs(co * 1E8 - 1459) < 1.0)
    print (f"Co, Spivey - matbal: {co}")

    co = Oil_Compressibility_Undersaturated(2993, 4000, 981, 41.5, .782, 184, \
         correlation = Correlation_Oil_Compressibility.Spivey_Valko_McCain, \
         application = Compressibility_Application.pta_modeling)
    assert(abs(co * 1E8 - 1541) < 1.0)
    print (f"Co, Spivey - pta: {co}")

def Oil_Compressibility_Saturated_Tests():
    co = Oil_Compressibility_Saturated(1591.7, 614.7, 185, 25.5, 0.76464, 0.7173, 1.12712, 0.00572235, 272)
    #assert(abs(co * 1E8 - 1668) < 1.0)
    print (f"Co: {co}")

def Oil_Density_Tests():
    r = Oil_Density(2993, 2993, 981, 41.5, 0.782, 0.834, 184)
    # assert(abs(co * 1E8 - 1668) < 1.0)
    print (f"Density, Spivey, pb: {r}")

    r = Oil_Density(2993, 1200, 441, 41.5, 0.782, 0.834, 184)
    #assert(abs(co * 1E8 - 1668) < 1.0)
    print (f"Density, Spivey, 1200: {r}")

    r = Oil_Density(2993, 4000, 981, 41.5, 0.782, 0.834, 184)
    #assert(abs(co * 1E8 - 1668) < 1.0)
    print (f"Density, Spivey, 4000: {r}")

def Oil_FVF_Tests():
    bo = Oil_FVF(2993, 2993, 981, 41.5, 0.834, correlation = Correlation_Oil_FVF.McCain, Density = 39.63)
    assert(abs(bo - 1.566) < 0.01)
    print (f"Bo: {bo}")

    bo = Oil_FVF(2993, 1200, 441, 41.5, 0.834, correlation = Correlation_Oil_FVF.McCain, Density = 43.68)
    assert(abs(bo - 1.282) < 0.01)
    print (f"Bo: {bo}")

    bo = Oil_FVF(2993, 4000, 981, 41.5, 0.834, correlation = Correlation_Oil_FVF.McCain, Density = 40.35)
    assert(abs(bo - 1.540) < 0.01)
    print (f"Bo: {bo}")

def Oil_Visocity_Tests():
    mo = Oil_Viscosity(2993, 2993, 981, 41.5, 184)
    assert(abs(mo - 0.335) < 0.001)
    print (f"Viscosity: {mo}")

def Oil_Reservoir_Gas_Gravity_Tests():
    ggres = Oil_Reservoir_Gas_Gravity(1200.0, 981.0, 41.5, 0.782, 184.0)
    #assert(abs(ggres - 0.781) < 0.001)
    print (f"Reservoir Gas Gravity: {ggres}")



# Bubble_Point_Tests()
# Solution_GOR_Tests()
# Oil_Compressibility_Tests()
#Oil_Density_Tests()
#Oil_FVF_Tests()
#Solution_GOR_Tests()
Oil_Compressibility_Saturated_Tests()
