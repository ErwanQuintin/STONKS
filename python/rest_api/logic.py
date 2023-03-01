'''
Created on 28 févr. 2023

@author: michel
'''
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

from core.STONKS_PreComputed_Position_alert import transient_alert, dec
from session_utils import SessionUtils

class ParamHolder:
    """
    structure carrying all source parameter in one shot
    """
    def __init__(self):
        self.ra = None
        self.dec = None
        self.pos_err = None
        self.flux = None
        self.flux_err = None
        self.band_flux = None
        self.band_fluxerr = None
        self.src_num = None
        self.pn_offax = None
        self.m1_offax = None
        self.m2_offax = None 

def process_one_observation(obsmli_path, queue): 
    print(f"Loading EPIC source list {obsmli_path}")
    try:
        raw_data = fits.open(obsmli_path, memmap=True)
    
        #Building Observation information using OBSMLI header
        dict_observation_metadata = {}
        obs_information = raw_data[0].header
        dict_observation_metadata["ObsID"] = str(obs_information['OBS_ID'])
        dict_observation_metadata["DateObs"] = obs_information['DATE-OBS']
        dict_observation_metadata["TargetName"] = obs_information['OBJECT']
        dict_observation_metadata["MJD"] = Time(dict_observation_metadata["DateObs"], format="isot").mjd
    
        sources_raw = raw_data[1].data
        sources_raw = Table(sources_raw)
        session_path = SessionUtils.get_session_path(dict_observation_metadata["ObsID"])
        tab_band_fluxes = [[list(line)] for line in sources_raw["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"]]
        tab_band_fluxerr = []
        nb_src = 0;
        nb_alerts = 0
        for line in sources_raw["ERR_EP_1_FLUX","ERR_EP_2_FLUX","ERR_EP_3_FLUX","ERR_EP_4_FLUX","ERR_EP_5_FLUX"]:
            tab_band_fluxerr.append([[list(line)], [list(line)]])
    
        for ra, dec, pos_err, flux, flux_err, band_flux, band_fluxerr, src_num, pn_offax, m1_offax, m2_offax in \
                zip(sources_raw["RA"],sources_raw["DEC"],sources_raw["RADEC_ERR"], sources_raw["EP_TOT_FLUX"],\
                        sources_raw["ERR_EP_TOT_FLUX"],tab_band_fluxes, tab_band_fluxerr, sources_raw["SRC_NUM"],
                    sources_raw["PN_OFFAX"],sources_raw["M1_OFFAX"],sources_raw["M2_OFFAX"]):
            param_holder = ParamHolder()
            param_holder.ra = ra 
            param_holder.dec = dec
            param_holder.pos_err = pos_err
            param_holder.flux = flux
            param_holder.flux_err = flux_err
            param_holder.band_flux = band_flux
            param_holder.band_fluxerr = band_fluxerr
            param_holder.src_num = src_num
            param_holder.pn_offax = pn_offax
            param_holder.m1_offax = m1_offax
            param_holder.m2_offax = m2_offax
            print (f"Processing source {src_num}")
            nb_alerts += process_one_source(param_holder, dict_observation_metadata) 
            nb_src += 1
            if nb_alerts  > 0:
                break
        if queue is not None:
            queue.put({"status": "succeed",
                       "session_name": dict_observation_metadata["ObsID"],
                       "nb_sources": str(nb_src),
                       "nb_alerts": str(nb_alerts),
                       })   
    except Exception as exp:
        print(exp)
        if queue is not None:
            queue.put({"status": "failed", "exception": f"{str(exp)}"})   
        
         
def process_one_source(param_holder, observation_metadata):
    tab_alerts=[]
    tab_dic_infos = []
    nb_alerts = 0
    dict_detection_info={}
    dict_detection_info["ObsID"]=observation_metadata["ObsID"]
    dict_detection_info["Date Obs"]=observation_metadata["DateObs"]
    dict_detection_info["Target Name"]=observation_metadata["TargetName"]
    dict_detection_info['SRCNUM'] = str(param_holder.src_num)
    dict_detection_info['Off-axis Angles'] = f'PN: {param_holder.pn_offax:.1f}", M1: {param_holder.m1_offax:.1f}", M2: {param_holder.m2_offax:.1f}"'
    dict_detection_info['Source RA']=np.round(param_holder.ra, 4)
    dict_detection_info['Source Dec']=np.round(param_holder.dec, 4)
    dict_detection_info['Position Error']=f'{param_holder.pos_err:.2f}"'

    result_alert = transient_alert(1,
                                   param_holder.ra,
                                   param_holder.dec,
                                   param_holder.pos_err,
                                   param_holder.flux,
                                   param_holder.flux_err,
                                   param_holder.band_flux,
                                   param_holder.band_fluxerr,
                                   observation_metadata["MJD"],
                                   var_flag=False)
    if result_alert!=[]: #If there is a transient alert
        tab_alerts += result_alert
        tab_dic_infos.append(dict_detection_info)
        for ms, dict_det_info in zip(tab_alerts,tab_dic_infos):
            ms.save_lightcurve(dict_det_info)
            nb_alerts += 1
 
    return nb_alerts