import numpy as np

def hrus_centroid_distance(lats,lons):
    radius = 6371 # km

    ncells = len(lats)
    distance = np.zeros((ncells,ncells))
    for i, lat1, lon1 in zip(range(ncells),lats,lons):
     for j, lat2, lon2 in zip(range(ncells),lats,lons):
      if j>i:
       dlat = np.radians(lat2-lat1)
       dlon = np.radians(lon2-lon1)
       a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(lat1)) \
         * np.cos(np.radians(lat2)) * np.sin(dlon/2) * np.sin(dlon/2)
       c = 2 * np.math.atan2(np.sqrt(a), np.sqrt(1-a))
       distance[i,j] = radius * c
       distance[j,i] = radius * c

    return distance

def hrus_slope(elevation,dist):
    ncells = len(elevation)

    slope = np.zeros((ncells,ncells))
    distance = np.array(dist)*1000.0 # convert to m

    for i, h1 in zip(range(ncells),elevation):
     for j, h2 in zip(range(ncells),elevation):
      if j>i : 
        slope[i,j] = (h1-h2)/distance[i,j]
        slope[j,i] = -slope[i,j]

    return slope

def calc_calendar(self,ncells):
    grow = np.zeros((ncells,12),dtype=int)
    st_gscal = np.copy(self.st_gscal)-1
    en_gscal = np.copy(self.en_gscal)-1    

    leng = en_gscal-st_gscal+1
   
    m = leng > 0
    for i in np.where(m)[0]:
      grow[i,st_gscal[i]:en_gscal[i]+1] = 1

    m = leng < 0
    for i in np.where(m)[0]:
      grow[i,:] = 1
      grow[i,en_gscal[i]+1:st_gscal[i]] = 0
 
  
    return grow



'''
def Simple_Local_Abstraction(self,NOAH,TOPMODEL):

  # Groundwater Abstractions
  if self.hwu_gw_flag == True :
    local_supply_gw = np.copy(self.supply_gw)

    # Domestic
    if self.hwu_domest_flag  == True:
        m = (local_supply_gw >= self.deficit_domest) & (self.deficit_domest > 0.0)
        local_supply_gw[m] = (local_supply_gw-self.deficit_domest)[m]
        self.deficit_domest[m] = 0.0
        m = (self.deficit_domest > local_supply_gw) & (local_supply_gw > 0.0)
        self.deficit_domest[m] = (self.deficit_domest-local_supply_gw)[m]
        local_supply_gw[m] = 0.0

    # Agriculture
    if self.hwu_agric_flag  == True:
        m = (local_supply_gw >= self.deficit_agric) & (self.deficit_agric > 0.0)
        local_supply_gw[m] = (local_supply_gw-self.deficit_agric)[m]
        self.deficit_agric[m] = 0.0
        m = (self.deficit_agric > local_supply_gw) & (local_supply_gw > 0.0)
        self.deficit_agric[m] = (self.deficit_agric-local_supply_gw)[m]
        local_supply_gw[m] = 0.0

    # Industrial
    if self.hwu_indust_flag == True:
        m = (local_supply_gw >= self.deficit_indust) & (self.deficit_indust > 0.0)
        local_supply_gw[m] = (local_supply_gw-self.deficit_indust)[m]
        self.deficit_indust[m] = 0.0
        m = (self.deficit_indust > local_supply_gw) & (local_supply_gw > 0.0)
        self.deficit_indust[m] = (self.deficit_indust-local_supply_gw)[m]
        local_supply_gw[m] = 0.0

    # Abstract 
    m = ((self.supply_gw-local_supply_gw) > 0.0)
    NOAH.dzwt[m] = (NOAH.dzwt-(self.supply_gw-local_supply_gw)/TOPMODEL.area)[m]
    self.supply_gw = local_supply_gw

  # Surface Water Abstraction
  if self.hwu_sf_flag == True :
    local_supply_sf = np.copy(self.supply_sf)

    # Domestic
    if self.hwu_domest_flag  == True:
        m = (local_supply_sf >= self.deficit_domest) & (self.deficit_domest > 0.0)
        local_supply_sf[m] = (local_supply_sf-self.deficit_domest)[m]
        self.deficit_domest[m] = 0.0
        m = (self.deficit_domest > local_supply_sf) & (local_supply_sf > 0.0)
        self.deficit_domest[m] = (self.deficit_domest-local_supply_sf)[m]
        local_supply_sf[m] = 0.0

    # Agriculture
    if self.hwu_agric_flag  == True:
        m = (local_supply_sf >= self.deficit_agric) & (self.deficit_agric > 0.0)
        local_supply_sf[m] = (local_supply_sf-self.deficit_agric)[m]
        self.deficit_agric[m] = 0.0
        m = (self.deficit_agric > local_supply_sf) & (local_supply_sf > 0.0)
        self.deficit_agric[m] = (self.deficit_agric-local_supply_sf)[m]
        local_supply_sf[m] = 0.0

    # Industrial
    if self.hwu_indust_flag == True:
        m = (local_supply_sf >= self.deficit_indust) & (self.deficit_indust > 0.0)
        local_supply_sf[m] = (local_supply_sf-self.deficit_indust)[m]
        self.deficit_indust[m] = 0.0
        m = (self.deficit_indust > local_supply_sf) & (local_supply_sf > 0.0)
        self.deficit_indust[m] = (self.deficit_indust-local_supply_sf)[m]
        local_supply_sf[m] = 0.0

    # Abstract 
    m = ((self.supply_sf-local_supply_sf) > 0.0)
    NOAH.runsf[m] = (NOAH.runsf - (self.supply_sf-local_supply_sf)*1000./NOAH.dt/TOPMODEL.area)[m]
    self.supply_sf = local_supply_sf

  return
'''



