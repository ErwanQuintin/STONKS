'''
Created on 28 févr. 2023

@author: michel
'''
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

from core.STONKS_PreComputed_Position_alert import transient_alert

def process_one_observation(obsmli_path): 
    print(f"Loading EPIC source list {obsmli_path}")
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

    tab_band_fluxes = [[list(line)] for line in sources_raw["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"]]
    tab_band_fluxerr = []
    for line in sources_raw["ERR_EP_1_FLUX","ERR_EP_2_FLUX","ERR_EP_3_FLUX","ERR_EP_4_FLUX","ERR_EP_5_FLUX"]:
        tab_band_fluxerr.append([[list(line)], [list(line)]])

    for ra, dec, pos_err, flux, flux_err, band_flux, band_fluxerr, src_num, pn_offax, m1_offax, m2_offax in \
            zip(sources_raw["RA"],sources_raw["DEC"],sources_raw["RADEC_ERR"], sources_raw["EP_TOT_FLUX"],\
                    sources_raw["ERR_EP_TOT_FLUX"],tab_band_fluxes, tab_band_fluxerr, sources_raw["SRC_NUM"],
                sources_raw["PN_OFFAX"],sources_raw["M1_OFFAX"],sources_raw["M2_OFFAX"]):
        print (f"Processing source {src_num}")
        process_one_source(ra, dec, pos_err, flux, flux_err, band_flux, band_fluxerr, src_num, pn_offax, m1_offax, m2_offax, dict_observation_metadata)          
        
         
def process_one_source(ra, dec, pos_err, flux, flux_err, band_flux, band_fluxerr, src_num, pn_offax, m1_offax, m2_offax, observation_metadata):
    tab_alerts=[]
    tab_dic_infos = []

    dict_detection_info={}
    dict_detection_info["ObsID"]=observation_metadata["ObsID"]
    dict_detection_info["Date Obs"]=observation_metadata["DateObs"]
    dict_detection_info["Target Name"]=observation_metadata["TargetName"]
    dict_detection_info['SRCNUM'] = str(src_num)
    dict_detection_info['Off-axis Angles'] = f'PN: {pn_offax:.1f}", M1: {m1_offax:.1f}", M2: {m2_offax:.1f}"'
    dict_detection_info['Source RA']=np.round(ra, 4)
    dict_detection_info['Source Dec']=np.round(dec, 4)
    dict_detection_info['Position Error']=f'{pos_err:.2f}"'

    result_alert = transient_alert(1, ra, dec, pos_err, flux, flux_err, band_flux,
                                 band_fluxerr, observation_metadata["MJD"], var_flag=False)
    if result_alert!=[]: #If there is a transient alert
        tab_alerts += result_alert
        tab_dic_infos.append(dict_detection_info)
        for ms, dict_det_info in zip(tab_alerts,tab_dic_infos):
            ms.save_lightcurve(dict_det_info)
    else :
        print("=> No variability detected")