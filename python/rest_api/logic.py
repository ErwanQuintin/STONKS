'''
Created on 28 févr. 2023

@author: michel
'''
import traceback
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord

from core.STONKS_PreComputed_Position_alert import transient_alert
from session import Session

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

def process_one_observation(session, obsmli_path, queue): 
    print(f"Loading EPIC source list {obsmli_path}")
    try:
        
        raw_data = fits.open(obsmli_path, memmap=True)
    
        #Building Observation information using OBSMLI header
        dict_observation_metadata = {}
        obs_information = raw_data[0].header
        dict_observation_metadata["ObsID"] = str(obs_information['OBS_ID'])
        dict_observation_metadata["DateObs"] = obs_information['DATE-OBS']
        dict_observation_metadata["TargetName"] = obs_information['OBJECT'].replace("_","\_")
        dict_observation_metadata["MJD"] = Time(dict_observation_metadata["DateObs"], format="isot").mjd

        sources_raw = raw_data[1].data
        sources_raw = Table(sources_raw)
        indices_not_spurious = (((sources_raw["PN_DET_ML"]>10) | (np.isnan(sources_raw["PN_DET_ML"]))) &
                                        ((sources_raw["M1_DET_ML"]>10) | (np.isnan(sources_raw["M1_DET_ML"]))) &
                                        ((sources_raw["M2_DET_ML"]>10) | (np.isnan(sources_raw["M2_DET_ML"]))))
        sources_raw = sources_raw[(sources_raw["EP_EXT_ML"]<6) & indices_not_spurious]
        tab_band_fluxes = [[list(line)] for line in sources_raw["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"]]
        tab_band_fluxerr = []
        nb_src = 0
        nb_alerts = 0
        for line in sources_raw["ERR_EP_1_FLUX","ERR_EP_2_FLUX","ERR_EP_3_FLUX","ERR_EP_4_FLUX","ERR_EP_5_FLUX"]:
            tab_band_fluxerr.append([[list(line)], [list(line)]])
    
        for ra, dec, pos_err, flux, flux_err, band_flux, band_fluxerr, src_num, pn_offax, m1_offax, m2_offax, pn_detml,m1_detml, m2_detml, ep_detml\
                in zip(sources_raw["RA"],sources_raw["DEC"],sources_raw["RADEC_ERR"], sources_raw["EP_TOT_FLUX"],\
                        sources_raw["ERR_EP_TOT_FLUX"],tab_band_fluxes, tab_band_fluxerr, sources_raw["SRC_NUM"],
                    sources_raw["PN_OFFAX"],sources_raw["M1_OFFAX"],sources_raw["M2_OFFAX"], sources_raw["PN_DET_ML"],
                       sources_raw["M1_DET_ML"],sources_raw["M2_DET_ML"],sources_raw["EP_DET_ML"]):
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
            param_holder.pn_detml = pn_detml
            param_holder.m1_detml = m1_detml
            param_holder.m2_detml = m2_detml
            param_holder.ep_detml = ep_detml

            print (f"Processing source {src_num}")
            nb_alerts += process_one_source(param_holder, dict_observation_metadata, session) 
            nb_src += 1
        if queue is not None:
            queue.put({"status": "succeed",
                       "obsid": dict_observation_metadata["ObsID"],
                       "nb_sources": str(nb_src),
                       "nb_alerts": str(nb_alerts),
                       })   
    except Exception as exp:
        print(exp)
        traceback.print_exc()
        if queue is not None:
            queue.put({"status": "failed", "exception": f"{str(exp)}"})   
        
         
def process_one_source(param_holder, observation_metadata, session):
    tab_alerts=[]
    tab_dic_infos = []
    nb_alerts = 0
    dict_detection_info={}
    dict_detection_info["ObsID"]=observation_metadata["ObsID"]
    dict_detection_info["Date Obs"]=observation_metadata["DateObs"]
    dict_detection_info["Target Name"]=observation_metadata["TargetName"]
    dict_detection_info['SRCNUM'] = str(param_holder.src_num)
    dict_detection_info['Off-axis Angles'] = f"PN: {param_holder.pn_offax:.1f}', M1: {param_holder.m1_offax:.1f}', M2: {param_holder.m2_offax:.1f}'"
    dict_detection_info['Instruments DetML'] = f'PN: {param_holder.pn_detml:.1f}, M1: {param_holder.m1_detml:.1f}, M2: {param_holder.m2_detml:.1f}, EP: {param_holder.ep_detml:.1f}'
    c=SkyCoord(param_holder.ra*u.deg,param_holder.dec*u.deg).to_string('hmsdms', sep=":")
    dict_detection_info['Source RA']= f'{np.round(param_holder.ra, 4)}   /   {c.split(" ")[0]}'
    dict_detection_info['Source Dec']= f'{np.round(param_holder.dec, 4)}   /   {c.split(" ")[1]}'
    dict_detection_info['Position Error']=f'{param_holder.pos_err:.2f}"'

    result_alert, flag_alert, info_source = transient_alert(session,
                                   session.obsid,
                                   param_holder.ra,
                                   param_holder.dec,
                                   param_holder.pos_err,
                                   param_holder.flux,
                                   param_holder.flux_err,
                                   param_holder.band_flux,
                                   param_holder.band_fluxerr,
                                   observation_metadata["MJD"],
                                   var_flag=False,
                                   )
    if result_alert!=[]: #If there is a transient alert
        tab_alerts += result_alert
        tab_dic_infos.append(dict_detection_info)
        for ms, dict_det_info in zip(tab_alerts,tab_dic_infos):
            ms.save_lightcurve(dict_det_info, flag_alert)
            nb_alerts += 1
 
    return nb_alerts