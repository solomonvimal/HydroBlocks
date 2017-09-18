#import sys
#sys.path.append('Preprocessing')
#sys.path.append('Model')
import datetime
import HydroBloks as HB
import sys

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))

 return metadata

#Read in the metadata file
metadata_file = sys.argv[1]
metadata = Read_Metadata_File(metadata_file)
ncores = metadata['parallel_ncores']

#Define the dates
idate = datetime.datetime(metadata['startdate']['year'],
			   metadata['startdate']['month'],
			   metadata['startdate']['day'],0)
fdate = datetime.datetime(metadata['enddate']['year'],
			   metadata['enddate']['month'],
			   metadata['enddate']['day'],23)

#Define the info
hydrobloks_info = {
        'icatch':metadata['catchment_id'],
	'input_file':metadata['input_file'],
	'output_file':metadata['output_file'],
        #'soil_file':metadata['soil_file'],
        'workspace':metadata['workspace'],
	'surface_flow_flag':metadata['surface_flow_flag'],
	'subsurface_flow_flag':metadata['subsurface_flow_flag'],
    "hwu_flag":metadata['hwu_flag'],
      "hwu_sf_flag":metadata['hwu_sf_flag'],
      "hwu_gw_flag":metadata['hwu_gw_flag'],
      "hwu_agric_flag":metadata['hwu_agric_flag'],
      "hwu_domest_flag":metadata['hwu_domest_flag'],
      "hwu_indust_flag":metadata['hwu_indust_flag'],
      "hwu_lstock_flag":metadata['hwu_lstock_flag'],
	'dt':metadata['dt'],#seconds
	'dtt':metadata['dtt'],#seconds
	'dx':metadata['dx'],#meters
	'nsoil':metadata['nsoil'],
	'ncores':metadata['parallel_ncores'],
	'idate':idate,
	'fdate':fdate,
	'nclusters_nc':metadata['nhru_nc'],
	'nclusters_c':metadata['nhru_c'],
	'nclusters':metadata['nhru_nc'] + metadata['nhru_c'],
	'model_type':metadata['model_type'],
        'create_mask_flag':metadata['create_mask_flag'],
        'mkl_flag':metadata['mkl_flag']
	}

#Run the model
HB.Run_Model(hydrobloks_info)
