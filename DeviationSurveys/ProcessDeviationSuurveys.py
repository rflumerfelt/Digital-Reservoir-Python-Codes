import fme
import fmeobjects
import numpy as np
import numpy.ma as ma
import math


class FeatureProcessor(object):

    def __init__(self):
        global lstAttributes
        lstAttributes = []
        
    def input(self, feature):
        global WellId, SurveyId, lstAttributes, md
        
        WellId = feature.getAttribute('Well_ID')
        SurveyId = feature.getAttribute('Survey_ID')
        md = feature.getAttribute('Measured_Depth')
        incl = feature.getAttribute('Inclination')
        az = feature.getAttribute('Azimuth')
        lstAttributes.append((md, incl, az))
        
    def close(self):
        pass

    def process_group(self):
        global lstAttributes
        
        self.Get_Arrays()
        cl = md2-md1
        dogleg = np.arccos(np.minimum(1, np.sin(i1) * np.sin(i2) * np.cos(a2-a1) + np.cos(i1) * np.cos(i2)))
        dogleg_deg = dogleg * 180.0 / math.pi
        #dls = dogleg_deg * 100 / ma.masked_values(cl,0,copy=True)
        #rf = np.tan(dogleg/2.0) * 2 /  ma.masked_values(dogleg,0,copy=True)
        dls = dogleg_deg * 100 / ma.masked_values(cl,0)
        rf = np.tan(dogleg/2.0) * 2 /  ma.masked_values(dogleg,0)
        rf=rf.filled(1.0)
        dls=dls.filled(0.0)
        dns = (np.sin(i1)*np.cos(a1) + np.sin(i2)*np.cos(a2)) * rf * cl / 2.0
        dew = (np.sin(i1)*np.sin(a1) + np.sin(i2) * np.sin(a2)) * rf * cl / 2.0
        dtvd = (np.cos(i1) + np.cos(i2))* rf * cl / 2.0
        ns = np.cumsum(dns)
        ew = np.cumsum(dew)
        tvd = np.cumsum(dtvd)
        dist = np.sqrt(np.square(ns) + np.square(ew))

        #for idx in range(len(lstAttributes)):
        for idx in range(len(lstAttributes)):
            feature = fmeobjects.FMEFeature()
            feature.setAttribute('Well_ID', WellId)
            feature.setAttribute('Survey_ID', SurveyId)
            feature.setAttribute('Measured_Depth', lstAttributes[idx][0])
            feature.setAttribute('Inclination', incl[idx])
            feature.setAttribute('Azimuth', az[idx])
            feature.setAttribute('Survey_ID', SurveyId)
            if idx > 0:
                feature.setAttribute('N_S_Calc', ns[idx-1])
                feature.setAttribute('E_W_Calc', ew[idx-1])
                feature.setAttribute('TVD_Calc', tvd[idx-1])
                feature.setAttribute('Distance_Projection_Calc', dist[idx-1])
                feature.setAttribute('Dog_Leg_Severity_Calc', dls[idx-1])
            else:
                feature.setAttribute('N_S_Calc', 0.0)
                feature.setAttribute('E_W_Calc', 0.0)
                feature.setAttribute('TVD_Calc', 0.0)
                feature.setAttribute('Distance_Projection_Calc', 0.0)
                feature.setAttribute('Dog_Leg_Severity_Calc', 0.0)
            self.pyoutput(feature)
        lstAttributes = []
        
    def Get_Arrays(self):
        global incl, i1, i2, az, a1, a2, md, md1, md2
        dt = np.dtype([('md', float),('incl', float), ('az', float)])
        data = np.array(lstAttributes, dtype=dt)
        np.sort(data,axis=0,order='md')
        incl = data['incl']
        i1 = data['incl'][:-1] * math.pi / 180.0
        i2 = data['incl'][1:] * math.pi / 180.0
        az = data['az']
        a1 = data['az'][:-1] * math.pi / 180.0
        a2 = data['az'][1:] * math.pi / 180.0
        md = data['md']
        md1 = data['md'][:-1]
        md2 = data['md'][1:]