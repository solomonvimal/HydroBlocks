import numpy as np
import os
from datetime import datetime
from scipy.sparse import csc_matrix, csr_matrix
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/pyHWU/')
import cost_functions_allocation as cost_funcs
import matplotlib.pyplot as plt


class Human_Water_Use:

 def __init__(self,HB,info):

  NOAH = HB.noahmp
  TOPMODEL = HB.dtopmodel
     
  ncells = NOAH.ncells
  self.itime = 0

  self.hwu_flag = info['hwu_flag']

  # Water Supply Flags
  self.hwu_gw_flag = info['hwu_gw_flag']
  self.hwu_sf_flag = info['hwu_sf_flag']
  if self.hwu_flag == True:
   if not any([self.hwu_gw_flag,self.hwu_sf_flag]): 
    exit('Error: All water supply flags are off - at least one water source is needed')

  # Water Demand Flags
  self.hwu_agric_flag  = info['hwu_agric_flag']
  self.hwu_domest_flag = info['hwu_domest_flag']
  self.hwu_indust_flag = info['hwu_indust_flag']
  self.hwu_lstock_flag = info['hwu_lstock_flag']
  demands = [self.hwu_agric_flag,self.hwu_domest_flag,self.hwu_indust_flag,self.hwu_lstock_flag]
  if self.hwu_flag == True:
   if not any(demands):
     exit('Error: All water demands flags are off - at least one water demand is needed') 

  # Land Cover for Water Use -- Must be the same defined at Preprocessing.py
  self.wuse_index = {}
  self.water_use_land_cover = {
    'industrial': [13,16,20],
    'domestic':   [7,8,9,10,13,16,19,20],
    'livestock':  [7,8,9,10,16,19,20],
    'agriculture':[12,14],
    'surface_water':[11,17]}
 
  # Optimal Allocation
  self.hwu_optimal_allocation_flag = True

  # GROUNDWATER supply variables
  if self.hwu_gw_flag == True:
   # Well_depht
   self.well_depth = -3
   # Correct well_depht for shallow soils and add 25% buffer.
   self.well_depth = np.maximum(self.well_depth,NOAH.zsoil[:,-1]*0.75)
   self.supply_gw = np.zeros(ncells,dtype=np.float64)
   self.alloc_gw = np.zeros(ncells,dtype=np.float64)

   # Groundwater pumping where clay content < 50% OR hydraulic conductivity > 10^-6 
   self.mask_gw = np.array([ False if i<=10.0**-6 or j>=50. else True for i,j in zip(NOAH.satdk0,NOAH.clay_pct) ])
   # Maximum groundwater pumping distance (km)
   self.gw_dist_lim  = 5.0 #km
   # Maximum slope allowed for upstream pumping
   self.gw_slope_lim = -1./1000. # m/m


  # SURFACE WATER supply variables
  if self.hwu_sf_flag == True:
   self.supply_sf = np.zeros(ncells,dtype=np.float64)
   self.alloc_sf  = np.zeros(ncells,dtype=np.float64)

   # Surface Water availability where ti > 10 and slope < 0.01
   cond1 = ( TOPMODEL.sti >= 10.0 ) & ( TOPMODEL.beta < 0.01)
   # Surface Water in water Bodies and/or wetlands
   cond2 = np.array([ True if x in self.water_use_land_cover["surface_water"] else False for x in NOAH.vegtyp ])
   self.mask_sf    = np.array([ True if i==True or j==True else False for i,j in zip(cond1,cond2) ])
   # Maximum surface water abstractions distance (km)
   self.sf_dist_lim  = 10.0 # km
   # Maximum slope allowed for upstream abstraction
   self.sf_slope_lim = -1./1000. # m/m 


  # AGRICULTURE demand variables
  if self.hwu_agric_flag == True:
   self.deficit_agric = np.zeros(ncells,dtype=np.float64)
   self.demand_agric  = np.zeros(ncells,dtype=np.float64)
   self.alloc_agric   = np.zeros(ncells,dtype=np.float64)
   self.irrigation    = np.zeros(ncells,dtype=np.float64)
   self.mask_agric    = [True if x in self.water_use_land_cover['agriculture'] else False for x in NOAH.vegtyp ]
   self.wuse_index['a'] = len(self.wuse_index)

   # Irrigation
   # 1: non-paddy irrigation, 2: paddy irrigation, 0: others
   self.irrig_land = HB.input_fp.groups['parameters'].variables['irrig_land'][:]
   self.mask_irrig = ( self.mask_agric == True) & ( self.irrig_land > 0.0 )
   self.mask_irrig = np.copy(self.mask_agric)

   # Crop calendar
   self.st_gscal = np.asarray(HB.input_fp.groups['parameters'].variables['start_growing_season'][:],dtype=np.int)
   self.en_gscal = np.asarray(HB.input_fp.groups['parameters'].variables['end_growing_season'][:],dtype=np.int)
   self.gscal = cost_funcs.calc_calendar(self,ncells)
   m = np.where(np.invert(self.mask_agric))[0]
   self.gscal[m,:] = 0
   # Test for the case with demand but zero irrigation


  # INDUSTRIAL demand variables
  if self.hwu_indust_flag == True:
   self.demand_indust  = np.zeros(ncells,dtype=np.float64)
   self.deficit_indust = np.zeros(ncells,dtype=np.float64)
   self.mask_indust    = [True if x in self.water_use_land_cover['industrial'] else False for x in NOAH.vegtyp ]
   self.wuse_index['i'] = len(self.wuse_index)


  # DOMESTIC demand variables
  if self.hwu_domest_flag == True:
   self.demand_domest  = np.zeros(ncells,dtype=np.float64)
   self.deficit_domest = np.zeros(ncells,dtype=np.float64)
   self.mask_domest    = [True if x in self.water_use_land_cover['domestic'] else False for x in NOAH.vegtyp ]
   self.wuse_index['d'] = len(self.wuse_index)


  # LIVESTOCK demand variables
  if self.hwu_lstock_flag == True:
   self.demand_lstock  = np.zeros(ncells,dtype=np.float64)
   self.deficit_lstock = np.zeros(ncells,dtype=np.float64)
   self.mask_lstock    = [True if x in self.water_use_land_cover['livestock'] else False for x in NOAH.vegtyp ]
   self.wuse_index['l'] = len(self.wuse_index)

  self.nwuse_index = len(self.wuse_index)



 def initialize_allocation(self,HB):
  NOAH = HB.noahmp
  TOPMODEL = HB.dtopmodel
  
 
  if self.hwu_flag == True:
   ncells = NOAH.ncells
   self.dta = HB.dt
   self.ntt = int(HB.dt/HB.dtt)

   # HRU distances
   # OPT 1) Centroid Distances
   #self.ctrd_lats = info['input_fp'].groups['parameters'].variables['centroid_lats'][:]
   #self.ctrd_lons = info['input_fp'].groups['parameters'].variables['centroid_lons'][:]
   #self.hrus_distances = cost_funcs.hrus_centroid_distance(self.ctrd_lats,self.ctrd_lons)  # km
   # OPT 2) Minimum Distance between a HRU centroid and it's neighboors boundary 
   self.hru_min_dist = HB.input_fp.groups['parameters'].variables['hru_min_dist'][:]
   self.hrus_distances = self.hru_min_dist  # km
   print "HRU's distances - mean:%f and std:%f" % (np.mean(self.hrus_distances[self.hrus_distances>0.]), np.std(self.hrus_distances[self.hrus_distances>0.]))
   #plt.hist(self.hru_min_dist.flatten())
   #plt.show()

   # HRU relative distances - Review this
   self.hrus_rel_dist = np.copy(self.hrus_distances)
   for i in range(ncells): self.hrus_rel_dist[i,:] = self.hrus_distances[i,:]/np.sum(self.hrus_distances[i,:])

   # Calculate slope between HRUs
   self.hrus_slopes = cost_funcs.hrus_slope(TOPMODEL.dem,self.hrus_distances)  # m/m


   # SURFACE WATER RATIO AND COST
   if self.hwu_sf_flag == True:
    # Surface Water Ratio 
    self.ratio_sf = np.ones((self.nwuse_index,ncells,ncells))

    # Update Ratio for Demands
    if self.hwu_agric_flag == True:
      m = np.invert(self.mask_irrig)
      self.ratio_sf[self.wuse_index['a'],:,m] = 0.0
    if self.hwu_domest_flag == True:
      m = np.invert(self.mask_domest)
      self.ratio_sf[self.wuse_index['d'],:,m] = 0.0
    if self.hwu_indust_flag == True:
      m = np.invert(self.mask_indust)
      self.ratio_sf[self.wuse_index['i'],:,m] = 0.0
    if self.hwu_lstock_flag == True:
      m = np.invert(self.mask_lstock)
      self.ratio_sf[self.wuse_index['l'],:,m] = 0.0

    # Update Ratio for Supply
    m = np.invert(self.mask_sf);                   self.ratio_sf[:,m,:] = 0.0
    m = (self.hrus_distances > self.sf_dist_lim ); self.ratio_sf[:,m] = 0.0
    m = (self.hrus_slopes    < self.sf_slope_lim); self.ratio_sf[:,m] = 0.0
    # sum over the diff water sectors, then over the demand nodes
    total_ratio_sf = np.sum(np.sum(self.ratio_sf, axis=0),axis=1)
    # update mask to only simulate over the active supply nodes
    self.mask_sf = (self.mask_sf == True) & ( total_ratio_sf > 0.0)

    # Update Ratio to relative values
    #m = self.ratio_sf > 0.0
    #for n in ln: 
    # for i in li:
    #   self.ratio_sf[n,i,:] = 1.0-self.hrus_distances[i,:]/np.sum(self.hrus_distances[i,:])

    # Surface Water Cost
    self.cost_sf    = np.ones((self.nwuse_index,ncells,ncells))
    self.cost_sf[:] = self.hrus_distances


   # GROUNDWATER RATIO AND COST 
   if self.hwu_gw_flag == True:

    # Groundwater Ratio
    self.ratio_gw = np.ones((self.nwuse_index,ncells,ncells))

    # Update Ratio for Demands
    if self.hwu_agric_flag == True:
      m = np.invert(self.mask_irrig)
      self.ratio_gw[self.wuse_index['a'],:,m] = 0.0
    if self.hwu_domest_flag == True:
      m = np.invert(self.mask_domest)
      self.ratio_gw[self.wuse_index['d'],:,m] = 0.0
    if self.hwu_indust_flag == True:
      m = np.invert(self.mask_indust)
      self.ratio_gw[self.wuse_index['i'],:,m] = 0.0
    if self.hwu_lstock_flag == True:
      m = np.invert(self.mask_lstock)
      self.ratio_gw[self.wuse_index['l'],:,m] = 0.0

    # Update Ratio for Supply
    m = np.invert(self.mask_gw);                   self.ratio_gw[:,m,:] = 0.0
    m = (self.hrus_distances > self.gw_dist_lim ); self.ratio_gw[:,m] = 0.0
    m = (self.hrus_slopes    < self.gw_slope_lim); self.ratio_gw[:,m] = 0.0
    # sum over the diff water sectors, then over the demand nodes
    total_ratio_gw = np.sum(np.sum(self.ratio_gw, axis=0),axis=1)
    # update mask to only simulate over the active supply nodes
    self.mask_gw = (self.mask_gw == True) & ( total_ratio_gw > 0.0)

    # Update Ratio to relative values
    '''print self.hrus_distances
    print self.ratio_gw[0]
    for i in range(ncells):
      m = self.ratio_gw[0,i,:] > 0
      a = self.hrus_distances[i,:][m]
      b = np.sum(self.hrus_distances[i,:][m])
      c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
      print a,b, c
      self.ratio_gw[0,i,:][m] = 1.0-c
    '''

    # Groundwater Cost
    self.cost_gw     = np.ones((self.nwuse_index,ncells,ncells))
    self.cost_gw[:]  = self.hrus_distances#*(self.sf_dist_lim/self.gw_dist_lim)


   # INITIALIZE ALLOCATION MODEL
   if self.hwu_optimal_allocation_flag == True:

    from pywr.core import Model as wr_Model
    from pywr.core import Input as wr_Input
    from pywr.core import Output as wr_Output
    from pywr.core import Link as wr_Link

    # Define Model
    #self.ntwkm = wr_Model(start="2017-01-01", end="2017-01-01", timestep=1, solver='lpsolve')
    self.ntwkm = wr_Model(start="2017-01-01", end="2017-01-01", timestep=1, solver='glpk')

    # Define Network Supply and Demand Nodes    
    if self.hwu_agric_flag == True:
     a_nodes_names=[] #Agriculture
     m = self.mask_irrig
     for i in np.where(m)[0]:
      dem = wr_Output(self.ntwkm, name='a%i' % i, min_flow = 0.0, max_flow = 0.0, cost=-999)
      a_nodes_names.append('a%i' % i)

    if self.hwu_indust_flag == True:
     i_nodes_names=[] #Industrial
     m = self.mask_indust
     for i in np.where(m)[0]:
      dem = wr_Output(self.ntwkm, name='i%i' % i, min_flow = 0.0, max_flow = 0.0, cost=-999)
      i_nodes_names.append('i%i' % i)

    if self.hwu_domest_flag == True:
     d_nodes_names=[] #Domestic
     m = self.mask_domest
     for i in np.where(m)[0]:
      dem = wr_Output(self.ntwkm, name='d%i' % i, min_flow = 0.0, max_flow = 0.0, cost=-999)
      d_nodes_names.append('d%i' % i)

    if self.hwu_lstock_flag == True:
     l_nodes_names=[] #Livestock
     m = self.mask_lstock
     for i in np.where(m)[0]:
      dem = wr_Output(self.ntwkm, name='l%i' % i, min_flow = 0.0, max_flow = 0.0, cost=-999)
      l_nodes_names.append('l%i' % i)

    if self.hwu_sf_flag == True:
     s_nodes_names=[] #Suface Water
     m = self.mask_sf
     for i in np.where(m)[0]:
      surf = wr_Input( self.ntwkm, name='s%i' % i, min_flow = 0.0, max_flow = 0.0, cost=1)
      s_nodes_names.append('s%i' % i)

    if self.hwu_gw_flag == True:
     g_nodes_names=[] #Groundwater
     m = self.mask_gw
     for i in np.where(m)[0]:
      gw = wr_Input( self.ntwkm, name='g%i' % i, min_flow = 0.0, max_flow = 0.0, cost=1)
      g_nodes_names.append('g%i' % i)


    # Define Surface Nodes Connections
    if self.hwu_sf_flag == True:
     s_links_names = []
     m = self.ratio_sf > 0.0
     sf_valid_links = m.flatten()
     for n,i,j in zip(*np.where(m)):
      n = self.wuse_index.keys()[n]
      link_sf = wr_Link(self.ntwkm, name='s%i_%s%i' % (i,n,j), min_flow = 0.0, max_flow = 0.0)
      self.ntwkm.nodes['s%i'%i].connect(link_sf)
      link_sf.connect(self.ntwkm.nodes['%s%i'% (n,j)])
      s_links_names.append('s%i_%s%i' % (i,n,j))

    # Define Groundwater Nodes Connections
    if self.hwu_gw_flag == True:
     g_links_names = []
     m = self.ratio_gw > 0.0
     gw_valid_links = m.flatten()
     for n,i,j in zip(*np.where(m)):
      n = self.wuse_index.keys()[n]
      link_gw = wr_Link(self.ntwkm, name='g%i_%s%i' % (i,n,j), min_flow = 0.0, max_flow = 0.0)
      self.ntwkm.nodes['g%i'%i].connect(link_gw)
      link_gw.connect(self.ntwkm.nodes['%s%i'% (n,j)])
      g_links_names.append('g%i_%s%i' % (i,n,j))

    # Check if there is any source of water in the catchment
    self.valid_links = np.empty(ncells)
    self.valid_links[:] = False
    if self.hwu_sf_flag == True and any(sf_valid_links):
      self.valid_links = sf_valid_links
      if self.hwu_gw_flag == True and any(gw_valid_links):
        self.valid_links = np.concatenate([sf_valid_links,gw_valid_links])
    elif self.hwu_gw_flag == True or any(gw_valid_links):
        self.valid_links = gw_valid_links



    #Model Setup
    # Test for valid nodes and connections
    if any(self.valid_links):

      #self.ntwkm.check()
      self.ntwkm.setup()

      self.nodes_list = self.ntwkm.graph.nodes(data=True)
      #print self.nodes_list[:5] 

      n_nodes = len(self.nodes_list)
      #print n_nodes

      #print dir(def_list[0][0])
      #print nodes_list[0][0].max_flow    
      nodes_names = [ str(self.nodes_list[i][0].name) for i in range(n_nodes) ]
      nodes_types = [ str(self.nodes_list[i][0].__class__.__name__) for i in range(n_nodes) ]

      # Order the list - List containing the ordered position of each node and link
      # Demands, Supply, Groundwater, Link_i_j
      if self.hwu_agric_flag == True:
        self.a_nodes_position = []
        for i in a_nodes_names: self.a_nodes_position.append( nodes_names.index(i) )
      if self.hwu_domest_flag == True:
        self.d_nodes_position = []
        for i in d_nodes_names: self.d_nodes_position.append( nodes_names.index(i) )
      if self.hwu_indust_flag == True:
        self.i_nodes_position = []
        for i in i_nodes_names: self.i_nodes_position.append( nodes_names.index(i) )
      if self.hwu_lstock_flag == True:
        self.l_nodes_position = []
        for i in l_nodes_names: self.l_nodes_position.append( nodes_names.index(i) )

      if self.hwu_sf_flag == True:
        self.s_nodes_position = []
        for i in s_nodes_names: self.s_nodes_position.append( nodes_names.index(i) )
        self.s_links_position = []
        for i in s_links_names: self.s_links_position.append( nodes_names.index(i) )

      if self.hwu_gw_flag == True:
        self.g_nodes_position = []
        for i in g_nodes_names: self.g_nodes_position.append( nodes_names.index(i) )
        self.g_links_position = []
        for i in g_links_names: self.g_links_position.append( nodes_names.index(i) )

      #print len(s_nodes_names), len(g_nodes_names), len(s_links_names), len(g_links_names)
      #print n_nodes
      #print d_nodes_names[:10]
      #print self.d_nodes_position[:10]
      #print [ nodes_names[i] for i in self.d_nodes_position[:10] ]

      #self.ntwkm.graph.nodes[0][0].max_flow=9876
      #print self.ntwkm.graph.nodes[0][0].max_flow
      #print "Allocation Network Done!"





    return 



 def Calc_Human_Water_Demand_Supply(self,HB):
   NOAH = HB.noahmp
   TOPMODEL = HB.dtopmodel
   
   if self.hwu_agric_flag  == True:
    # Calculate Agricultural Demand
    self.demand_agric = self.Agriculture_Demand(HB)/self.dta #m/s
    # Convert from m/s to m3
    self.deficit_agric = np.copy(self.demand_agric)*self.dta*TOPMODEL.area

   if self.hwu_flag == True:
    # Convert from m/s to m3
    if self.hwu_indust_flag == True: self.deficit_indust = np.copy(self.demand_indust)*self.dta*TOPMODEL.area
    if self.hwu_domest_flag == True: self.deficit_domest = np.copy(self.demand_domest)*self.dta*TOPMODEL.area
    if self.hwu_lstock_flag == True: self.deficit_lstock = np.copy(self.demand_lstock)*self.dta*TOPMODEL.area

   if self.hwu_flag == True:
    # Calculate Supply [m]
    self.Calc_Water_Supply(HB)
    # Convert m to m3
    if self.hwu_sf_flag == True: self.supply_sf = self.supply_sf*TOPMODEL.area
    if self.hwu_gw_flag == True: self.supply_gw = self.supply_gw*TOPMODEL.area
  
    gvol = np.copy(self.supply_gw)[:,np.newaxis]+1.
    ncost = np.copy(self.hrus_distances+1)/gvol
    ncost[ncost == np.inf] = 10000000000.
    ncost[ncost == np.nan] = 10000000000.
    #test = self.hrus_distances[:,0]/self.supply_gw[0]
    #print test[:4]

    # Allocation
    if self.hwu_optimal_allocation_flag == True: 
     self.Optimal_Water_Allocation(NOAH,TOPMODEL)
   
   # Convert from m3 to m/s
   if self.hwu_agric_flag  == True: self.deficit_agric  = self.deficit_agric/self.dta/TOPMODEL.area

   if self.hwu_flag == True:
    # Convert from m3 to m/s  
    if self.hwu_indust_flag == True: self.deficit_indust = self.deficit_indust/self.dta/TOPMODEL.area
    if self.hwu_domest_flag == True: self.deficit_domest = self.deficit_domest/self.dta/TOPMODEL.area
    if self.hwu_lstock_flag == True: self.deficit_lstock = self.deficit_lstock/self.dta/TOPMODEL.area

    if self.hwu_agric_flag  == True: self.alloc_agric  = self.alloc_agric/TOPMODEL.area/self.ntt
    if self.hwu_gw_flag == True : self.alloc_gw = self.alloc_gw/TOPMODEL.area/self.ntt
    if self.hwu_sf_flag == True : self.alloc_sf = self.alloc_sf/TOPMODEL.area/self.ntt

    if self.hwu_sf_flag == True: self.supply_sf = self.supply_sf/TOPMODEL.area
    if self.hwu_gw_flag == True: self.supply_gw = self.supply_gw/TOPMODEL.area
   
   #print self.date, self.alloc_agric, self.alloc_gw
   return


 def Calc_Water_Supply(self,HB):
  NOAH = HB.noahmp
  TOPMODEL = HB.dtopmodel
  
  if self.hwu_flag == True:

   # Groundwater Supply
   if self.hwu_gw_flag == True :
    self.supply_gw = (NOAH.zwt-self.well_depth)*(NOAH.smcmax-NOAH.wltsmc0) # water height above well [m]
    m = (self.supply_gw < 0)
    self.supply_gw[m] = 0.0
    # Future: Include corrections for environmental flows
    self.supply_gw = self.supply_gw*0.8

   # Surface Water Supply
   if self.hwu_sf_flag == True :
    self.supply_sf = NOAH.runsf*self.dta/1000. # from mm/s to m
    # Future Include correction for environmetal flows.
    self.supply_sf = self.supply_sf*0.8

  return


 def Water_Supply_Abstraction(self, HB):
  NOAH = HB.noahmp
  TOPMODEL = HB.dtopmodel

  if self.hwu_flag == True:

   # Abstract from Surface
   if self.hwu_sf_flag == True:
    m = (self.alloc_sf > 0.0)
    NOAH.runsf[m] = (NOAH.runsf - (self.alloc_sf*1000/NOAH.dt))[m] # mm/s 

   # Abstract from Grondwater
   if self.hwu_gw_flag == True:
    m = (self.alloc_gw > 0.0)
    NOAH.dzwt[m] = (NOAH.dzwt-(self.alloc_gw))[m]  # m

 def Human_Water_Irrigation(self,HB):
  NOAH = HB.noahmp
  TOPMODEL = HB.dtopmodel
  
  if self.hwu_flag == True:
   if self.hwu_agric_flag  == True:
    # Add as irrigation the amount of water that was allocated
    self.irrigation[:] = self.alloc_agric*1000/NOAH.dt 
    #m = (self.irrigation > 0.0) & (self.mask_irrig == True)
    m = (self.mask_irrig == True)
    NOAH.prcp[m] = (NOAH.prcp + self.irrigation)[m]
    # Include irrigation efficiency later on

  return


 def Optimal_Water_Allocation(self,NOAH,TOPMODEL):  

    ncells = NOAH.ncells
   
    if any(self.valid_links.flatten()):

     p1 = datetime.now()
     #print 'P1', str(p1)
     # Update Node
     if self.hwu_agric_flag  == True:
      a_deficit = np.copy(self.deficit_agric)
      for p,d in zip(self.a_nodes_position,a_deficit[self.mask_irrig]):
        self.nodes_list[p][0].max_flow = d
     
     if self.hwu_indust_flag  == True:
      i_deficit = np.copy(self.deficit_indust)
      for p,d in zip(self.i_nodes_position,i_deficit[self.mask_indust]):
        self.nodes_list[p][0].max_flow = d
 
     if self.hwu_domest_flag  == True:
      d_deficit = np.copy(self.deficit_domest)
      for p,d in zip(self.d_nodes_position,d_deficit[self.mask_domest]):
        self.nodes_list[p][0].max_flow = d
  
     if self.hwu_lstock_flag  == True:
      l_deficit = np.copy(self.deficit_lstock)
      for p,d in zip(self.l_nodes_position,l_deficit[self.mask_lstock]):
        self.nodes_list[p][0].max_flow = d

     if self.hwu_sf_flag  == True:
      sf_supply = np.copy(self.supply_sf)
      for p,s in zip(self.s_nodes_position,sf_supply[self.mask_sf]):
        self.nodes_list[p][0].max_flow = s

     if self.hwu_gw_flag  == True:
      gw_supply = np.copy(self.supply_gw)
      for p,s in zip(self.g_nodes_position,gw_supply[self.mask_gw]):
        self.nodes_list[p][0].max_flow = s
      

     p2 = datetime.now()
     #print 'P2', str(p2-p1), str(p2-p1)
     # Update Links
 
     # Surface Water
     if self.hwu_sf_flag  == True:
      m = self.ratio_sf > 0.0
      cost  = self.cost_sf[m]
      max_flow = np.transpose(np.transpose(self.ratio_sf,(0, 2, 1))*sf_supply,(0,2,1))[m]
      for p,f,c in zip(self.s_links_position, max_flow, cost):
       self.nodes_list[p][0].max_flow = f
       self.nodes_list[p][0].cost = c  
    
     # Groundwater
     if self.hwu_gw_flag  == True:
      m = self.ratio_gw > 0
      cost  = self.cost_gw[m]
      max_flow =  np.transpose(np.transpose(self.ratio_gw,(0, 2, 1))*gw_supply,(0,2,1))[m]
      for p,f,c in zip(self.g_links_position, max_flow, cost):
       self.nodes_list[p][0].max_flow = f
       self.nodes_list[p][0].cost = c

     # Optimize the network
     total_deficit = 0.0; total_sf_supply= 0.0; total_gw_supply=0.0
     if self.hwu_agric_flag == True:  total_deficit +=np.sum(a_deficit[self.mask_irrig])
     if self.hwu_domest_flag == True: total_deficit +=np.sum(d_deficit[self.mask_domest])
     if self.hwu_indust_flag == True: total_deficit +=np.sum(i_deficit[self.mask_indust])
     if self.hwu_lstock_flag == True: total_deficit +=np.sum(l_deficit[self.mask_lstock])
     if self.hwu_sf_flag == True: total_sf_supply = np.sum(sf_supply[self.mask_sf])
     if self.hwu_gw_flag == True: total_gw_supply = np.sum(gw_supply[self.mask_gw])

     if total_deficit != 0.0 and ( total_sf_supply != 0.0 or total_gw_supply != 0.0 ): 
    
        p3 = datetime.now()
        #print 'P3', str(p3-p1), str(p3-p2)
        self.ntwkm.run()
  
        p4 = datetime.now()
        #print 'P4', str(p4-p1), str(p4-p3)
        
        delta = np.zeros(ncells)     
        
        if self.hwu_agric_flag == True:
         for p,h in zip(self.a_nodes_position,np.where(self.mask_irrig)[0]):
           delta[h] = self.nodes_list[p][0].flow[0]
         #print " " 
         #print a_deficit, delta
         #print " " 
         a_deficit = a_deficit - delta
         
        if self.hwu_indust_flag == True:
         delta[:] = 0     
         for p,h in zip(self.i_nodes_position,np.where(self.mask_indust)[0]):
           delta[h] = self.nodes_list[p][0].flow[0]
         i_deficit = i_deficit - delta

        if self.hwu_domest_flag == True:
         delta[:] = 0     
         for p,h in zip(self.d_nodes_position,np.where(self.mask_domest)[0]):
           delta[h] = self.nodes_list[p][0].flow[0]
         d_deficit = d_deficit - delta

        if self.hwu_lstock_flag == True:
         delta[:] = 0     
         for p,h in zip(self.l_nodes_position,np.where(self.mask_lstock)[0]):
           delta[h] = self.nodes_list[p][0].flow[0]
         l_deficit = l_deficit - delta
        
        if self.hwu_sf_flag == True:
         delta[:] = 0
         for p,h in zip(self.s_nodes_position,np.where(self.mask_sf)[0]):
           delta[h] = self.nodes_list[p][0].flow[0]
         sf_supply = sf_supply - delta
    
        if self.hwu_gw_flag == True:
         delta[:] = 0
         for p,h in zip(self.g_nodes_position,np.where(self.mask_gw)[0]):
           delta[h] = self.nodes_list[p][0].flow[0]
         gw_supply = gw_supply - delta

        p5 = datetime.now()
        #print 'P5', str(p5-p1), str(p5-p4)


     # Update Allocated water
     if self.hwu_agric_flag == True: self.alloc_agric = np.copy(self.deficit_agric) -a_deficit
     if self.hwu_sf_flag == True: self.alloc_sf = np.copy(self.supply_sf) -sf_supply
     if self.hwu_gw_flag == True: self.alloc_gw = np.copy(self.supply_gw) -gw_supply
     
     # Update Deficits 
     if self.hwu_agric_flag == True:  self.deficit_agric  = a_deficit
     if self.hwu_domest_flag == True: self.deficit_domest = d_deficit
     if self.hwu_indust_flag == True: self.deficit_indust = i_deficit
     if self.hwu_lstock_flag == True: self.deficit_lstock = l_deficit
    
     if self.hwu_sf_flag == True: self.supply_sf = sf_supply
     if self.hwu_gw_flag == True: self.supply_gw = gw_supply
     
     return


 def Agriculture_Demand(self,HB):
  NOAH = HB.noahmp
  TOPMODEL = HB.dtopmodel
  HWU = HB.hwu
  
  demand_agric = self.Calculate_Irrigation_Deficit(NOAH) # m
 
  # Limit Demand for the Crop Calendar
  mnt = HB.idate.month-1
  m = HWU.gscal[:,mnt] == 0
  demand_agric[m] = 0.0  

  return demand_agric


 def Calculate_Irrigation_Deficit(self,NOAH):
  """
  Calculates the irrigation deficit based on how much water [m] is necessary to 
  soil moisture reach field capacity
  """
  ncells = NOAH.ncells
  nroot_zone_depht = NOAH.root_depth
  smc = NOAH.sh2o
  sldpth = NOAH.sldpth
  smcref = NOAH.smcref # Field Capacity
  wltsmc0 = NOAH.wltsmc0 # Wilting point
  smcmax = NOAH.smcmax # Saturated 
  
  self.demand_agric[:] = 0.0

  # Irrigation Demand for general crop areas
  m = self.mask_agric
  if any(m):
   dsm = smcref[:,np.newaxis] - smc
   dsm [dsm < 0] = 0.0
   soil_demand = np.array([ np.sum(dsm[i,:nroot_zone_depht[i]]*sldpth[i,:nroot_zone_depht[i]]) for i in range(ncells) ])
   self.demand_agric[m] = soil_demand[m]
  
  # Irrigation Demand for paddy  crop areas
  m = ( self.irrig_land == 2.0 ) & (self.mask_agric == True)
  if any(m):
   dsm = smcmax[:,np.newaxis] - smc
   dsm[dsm < 0] = 0.0
   soil_demand = np.array([ np.sum(dsm[i,:nroot_zone_depht[i]]*sldpth[i,:nroot_zone_depht[i]]) for i in range(ncells) ])
   self.demand_agric[m] = soil_demand[m]  
  

  # Check if irrigation is triggered or not - Ozdodan et. al (2010)
  #MAi = (smc - wltsmc0[:,np.newaxis])/(sldpth-wltsmc0[:,np.newaxis])
  #MA = np.array([ np.sum(MAi[i,:nroot_zone_depht[i]]*sldpth[i,:nroot_zone_depht[i]])/np.sum(sldpth[i,:nroot_zone_depht[i]]) for i in range(ncells) ])
  #self.mask_trigger = (MA > 0.5*smcref) & (self.mask_irrig == True)
  
   

  return self.demand_agric


    

