import numpy as np
import os

class Human_Water_Use:

 def __init__(self,NOAH,ncells):
  self.ncells = ncells
  self.itime = 0
  self.well_depth = -3.0
  self.irrig_gw_supply = np.zeros(ncells,dtype=np.float64)
  self.irrig_sf_supply = np.zeros(ncells,dtype=np.float64)
  self.irrig_demand = np.zeros(ncells,dtype=np.float64)
  self.hwu_gw_flag = False
  self.hwu_sf_flag = True

  
 def Calculate_Irrigation_Deficit(self,smcref,sldpth,smc,nroot_zone_depht,ncells):
  """
  Calculates the irrigation deficit based on how much water [m] is necessary to 
  soil moisture reach field capacity
  """
  irrig_demand = np.zeros(ncells,dtype=np.float64)
  dsm = smcref[:,np.newaxis] - smc
  dsm [dsm < 0] = 0.0
  irrig_demand = [ np.sum(dsm[i,:nroot_zone_depht[i]]*sldpth[i,:nroot_zone_depht[i]])  for i in range(ncells) ]
  return irrig_demand


 def Irrigation_Demand(self,NOAH,HB):
  #Calculate irrigation demand
  self.irrig_demand[:] = self.Calculate_Irrigation_Deficit(NOAH.smcref,NOAH.sldpth,NOAH.sh2o,NOAH.root_depth,NOAH.ncells)
  return


 def Human_Water_Demand(self,NOAH,TOPMODEL,HB):
  # Irrigation Demand
  self.Irrigation_Demand(NOAH,HB)

  # Domestic Demand

  # Industrial Demand
  
  return
   

 def Human_Water_Abstraction(self,NOAH,TOPMODEL,HB):

  self.irrig_gw_supply[:] = 0.0
  self.irrig_sf_supply[:] = 0.0 

  # Groundwater avaialable supply and abstraction
  if self.hwu_gw_flag == True :
    # GW: check how much groundwater you have available
    water_height_above_well = (NOAH.zwt-self.well_depth)*NOAH.smcmax
    m = (self.irrig_demand > 0.0) & (water_height_above_well > 0.0)
    self.irrig_gw_supply[m] = self.irrig_demand[m]
    m = (self.irrig_gw_supply > water_height_above_well) & (water_height_above_well > 0.0)
    self.irrig_gw_supply[m] = water_height_above_well[m]

    # GW: Abstract water from available groundwater
    m = self.irrig_gw_supply > 0.0
    NOAH.dzwt[m] = (NOAH.dzwt-self.irrig_gw_supply)[m]

  # Check how much of the initial demand was not yet met by groundwater
  self.irrig_demand_defict = self.irrig_demand - self.irrig_gw_supply
  
  # Surface available supply and abstraction
  if self.hwu_sf_flag == True :
    # SF: check how much surface water you have available
    runoff_supply = NOAH.runsf*NOAH.dt/1000. # from mm/s to m
    m = (runoff_supply > 0.0) & (self.irrig_demand_defict > 0)
    self.irrig_sf_supply[m] = self.irrig_demand_defict[m]
    m = (self.irrig_demand_defict > runoff_supply) & (runoff_supply > 0.0)
    self.irrig_sf_supply[m] = runoff_supply[m] 

    # Abstract water from local available surface water
    m = self.irrig_sf_supply > 0.0
    NOAH.runsf[m] = (NOAH.runsf - self.irrig_sf_supply*1000./NOAH.dt)[m]

  # Abstract water from regional available surface water

  return


 def Human_Water_Irrigation(self,NOAH,TOPMODEL,HB):
  # Add as irrigation the amount of water that was allocated to supply the demand 
  if self.hwu_gw_flag == True :
    m = self.irrig_gw_supply > 0
    NOAH.prcp[m] = (NOAH.prcp+self.irrig_gw_supply*1000./NOAH.dt)[m]
  
  if self.hwu_sf_flag == True :
    m = self.irrig_sf_supply > 0
    NOAH.prcp[m] = (NOAH.prcp+self.irrig_sf_supply*1000./NOAH.dt)[m]

  return


