import fme
import fmeobjects
import numpy as np
import math

class FeatureProcessor(object):

    def __init__(self):
        params = self.getParameters(['MAX_SPACING', 'MAX_HEEL_TOE_DISTANCE','GROUP_BY','X1','Y1','X2','Y2', 'TVD'])
        self.max_spacing = float(params['MAX_SPACING'])
        self.max_heel_toe_dist = float(params['MAX_HEEL_TOE_DISTANCE'])
        self.well_id_name = params['GROUP_BY']
        self.x1_name = str(params['X1'])
        self.y1_name = str(params['Y1'])
        self.x2_name = str(params['X2'])
        self.y2_name = str(params['Y2'])
        self.tvd_name = str(params['TVD'])
        self.id=[]
        self.x1=[]
        self.y1=[]
        self.x2=[]
        self.y2=[]
        self.tvd=[]

    def input(self, feature):
        id = feature.getAttribute(self.well_id_name)
        print(id)
        x1 = feature.getAttribute(self.x1_name)
        y1 = feature.getAttribute(self.y1_name)
        x2 = feature.getAttribute(self.x2_name)
        y2 = feature.getAttribute(self.y2_name)
        tvd = feature.getAttribute(self.tvd_name)
        self.id.append(id)
        self.x1.append(float(x1))
        self.y1.append(float(y1))
        self.x2.append(float(x2))
        self.y2.append(float(y2))
        if tvd == "":
           tvd = 0
        self.tvd.append(float(tvd))

    def close(self):
         # Initialize variables
         max_heel_toe_dist = self.max_heel_toe_dist
         max_spacing = self.max_spacing
         id = np.array(self.id)
         x1 = np.array(self.x1)
         y1 = np.array(self.y1)
         x2 = np.array(self.x2)
         y2 = np.array(self.y2)
         tvdo = np.array(self.tvd)
         az = self.azimuth(x1,y1,x2,y2) #angle relative to due north
         n = len(id) #number of wells
         for idx in range(n):
             #rotate well(w) and offsets(o) so that the target well is due north south
             sin_az = math.sin(az[idx])
             cos_az = math.cos(az[idx])
             x1o = x1 * cos_az - y1 * sin_az
             y1o = x1 * sin_az + y1 * cos_az
             x2o = x2 * cos_az - y2 * sin_az
             y2o = x2 * sin_az + y2 * cos_az
             #get points at ymin and ymax values for well and offsets
             tvdo = np.array(self.tvd)
             yw_min = min(y1o[idx], y2o[idx])
             yw_max = max(y1o[idx], y2o[idx])
             xw = x1o[idx]
             tvdw = tvdo[idx]
             yo_min = np.copy(y1o)
             yo_max = np.copy(y1o)
             xo_min = np.copy(x1o)
             xo_max = np.copy(x1o)
             yo_min[y2o < y1o] = y2o[y2o < y1o]
             yo_max[y2o > y1o] = y2o[y2o > y1o]
             xo_min[y2o < y1o] = x2o[y2o < y1o]  # note this is the x corresponding to the min y, not the minimum x
             xo_max[y2o > y1o] = x2o[y2o > y1o]  # note this is the x corresponding to the max y, not the maximum x
             #Calculate heel-toe offset and filter
             ht1 = yo_min - yw_max
             ht2 = yw_min - yo_max
             fltr = np.logical_not(np.logical_or(ht1 > max_heel_toe_dist, ht2 > max_heel_toe_dist))
             heel_toe_dist = np.zeros(len(fltr[fltr]))
             ht1 = ht1[fltr]
             heel_toe_dist[ht1 > 0] = ht1[ht1 > 0]
             ht2 = ht2[fltr]
             heel_toe_dist[ht2 > 0] = ht2[ht2 > 0]
             #adjust offset endpoints beyond lateral
             idf = np.copy(id[fltr])
             xo_min = xo_min[fltr]
             xo_max = xo_max[fltr]
             yo_min = yo_min[fltr]
             yo_max = yo_max[fltr]
             yo_min_adj = np.copy(yo_min)
             yo_min_adj[yo_min < yw_min] = yw_min
             xo_min_adj = np.copy(xo_min)
             xo_min_adj[np.logical_and(yo_min < yw_min, np.absolute(yo_max-yo_min) > 0.0001) ] = (xo_min + (xo_max - xo_min) / (yo_max - yo_min)*(yw_min - yo_min))[yo_min < yw_min]
             yo_max_adj = np.copy(yo_max)
             yo_max_adj[yo_max > yw_max] = yw_max
             xo_max_adj = np.copy(xo_max)
             xo_max_adj[yo_max > yw_max] = (xo_min + (xo_max - xo_min) / (yo_max - yo_min)*(yw_max-yo_min))[yo_max > yw_max]
             spacing = (xo_max_adj + xo_min_adj) / 2 - xw
             is_left = np.array([True] * len(idf))
             is_left[spacing > 0] = False
             spacing = np.absolute(spacing)
             pct_overlap = (yo_max_adj - yo_min_adj) / (yw_max-yw_min) * 100
             tvdo = tvdo[fltr]
             fltr_spcg = spacing < max_spacing
             ido = idf[fltr_spcg]
             spacing = spacing[fltr_spcg]
             is_left = is_left[fltr_spcg]
             heel_toe_dist = heel_toe_dist[fltr_spcg]
             pct_overlap = pct_overlap[fltr_spcg]
             pct_overlap[heel_toe_dist > 0] = 0
             tvdo = tvdo[fltr_spcg]
             for j in range(len(ido)):
                 if is_left[j]:
                     lr = 'Left'
                 else:
                     lr = 'Right'
                 feature = fmeobjects.FMEFeature()
                 feature.setAttribute('Well_ID',int(id[idx]))
                 feature.setAttribute('Offset_Well_ID',int(ido[j]))
                 feature.setAttribute('Well_Spacing',float(spacing[j]))
                 feature.setAttribute('Left_Right',lr)
                 feature.setAttribute('Pct_Overlap',float(pct_overlap[j]))
                 feature.setAttribute('Heel_Toe_Dist',float(heel_toe_dist[j]))
                 if tvdo[j]!=0 and tvdw!=0:
                     feature.setAttribute('TVD_Difference',float(tvdo[j]-tvdw))
                 self.pyoutput(feature)

    def process_group(self):
         pass
         
    def azimuth(self, x1, y1, x2, y2):
         az = np.empty(len(x1))
         due_north = np.absolute(x1 - x2) < 0.0001
         not_due_north = np.logical_not(due_north)
         az[due_north] = 0
         az[np.logical_not(due_north)] = np.arctan((x1 - x2) / (y1 - y2))
         return az



    def getParameters(self,ParameterNames):
          ret = dict()
          for macroName in FME_MacroValues:
              for name in ParameterNames:
                  if name in str(macroName):
                      ret[name] = FME_MacroValues[str(macroName)]
                      break
          return ret