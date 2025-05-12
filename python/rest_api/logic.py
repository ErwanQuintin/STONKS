'''
Created on 28 févr. 2023

@author: michel
'''
import traceback
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
import json

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
        self.off_target_angle = None
        self.pn_detml = None
        self.m1_detml = None
        self.m2_detml = None
        self.ep_detml = None
        self.choice_pi = None
        self.publishable = None
        self.pointing_type = None


def process_one_observation(session, queue):
    print(f"Loading EPIC source list {session.obsmli_path}")
    try:
        if session.image_path is not None:
            #Loads the EPIC image if it was uploaded
            raw_data = fits.open(session.image_path, memmap=True)
            image_wcs = wcs.WCS(raw_data[0].header)
            image_data = raw_data[0].data
            ra_target, dec_target = raw_data[0].header['RA_OBJ'], raw_data[0].header['DEC_OBJ']
            target_position = SkyCoord(ra=ra_target, dec=dec_target, unit='deg')
            flag_pointing_type='target'
        else:
            image_data=None
            image_wcs =None

        #Loads the data from the sources
        raw_data = fits.open(session.obsmli_path, memmap=True)
        if 'VARALERT' in raw_data[0].header.keys():
            choice_PI = raw_data[0].header['VARALERT']
        else:
            print('keyword for PI choice not present - assume no publishable alert')
            choice_PI=3
        if session.image_path is None:
            #If no image was given, we use the pointing position as a proxy to the target position, NOT IDEAL
            ra_target, dec_target = raw_data[0].header['RA_PNT'], raw_data[0].header['DEC_PNT']
            target_position = SkyCoord(ra=ra_target, dec=dec_target, unit='deg')
            flag_pointing_type='pointing'

        min_off_axis_angle = 2 #Minimum accepted off-axis angle in arcmin, to reject the source
        min_det_ml = 10 #Minimum accepted detection likelihood

        #Building Observation information using OBSMLI header
        dict_observation_metadata = {}
        obs_information = raw_data[0].header
        dict_observation_metadata["ObsID"] = str(obs_information['OBS_ID'])
        dict_observation_metadata["DateObs"] = obs_information['DATE-OBS']
        dict_observation_metadata["TargetName"] = obs_information['OBJECT'].replace("_","\_")
        dict_observation_metadata["MJD"] = Time(dict_observation_metadata["DateObs"], format="isot").mjd
        dict_observation_metadata["ExpTime"] = np.max((int(obs_information['EXPOS_PN']),
                                                       int(obs_information['EXPOS_M1']),
                                                       int(obs_information['EXPOS_M2'])))


        sources_raw = raw_data[1].data
        sources_raw = Table(sources_raw)
        indices_not_spurious = (((sources_raw["PN_DET_ML"]>min_det_ml) | (np.isnan(sources_raw["PN_DET_ML"]))) &
                                        ((sources_raw["M1_DET_ML"]>min_det_ml) | (np.isnan(sources_raw["M1_DET_ML"]))) &
                                        ((sources_raw["M2_DET_ML"]>min_det_ml) | (np.isnan(sources_raw["M2_DET_ML"]))))
        #indices_not_target = (((sources_raw["PN_OFFAX"] > min_off_axis_angle) | (np.isnan(sources_raw["PN_OFFAX"]))) &
        #                        ((sources_raw["M1_OFFAX"] > min_off_axis_angle) | (np.isnan(sources_raw["M1_OFFAX"]))) &
        #                        ((sources_raw["M2_OFFAX"] > min_off_axis_angle) | (np.isnan(sources_raw["M2_OFFAX"]))))
        sources_positions = SkyCoord(ra=sources_raw['RA'], dec=sources_raw['DEC'],unit='deg')
        off_target_angles = target_position.separation(sources_positions).to(u.arcmin).value
        indices_not_target = off_target_angles>min_off_axis_angle
        indices_not_target=indices_not_target[(sources_raw["EP_EXT_ML"]<6) & indices_not_spurious]

        #Choice_pi is going to be one of three (1,2,3): everything publishable, just serendipitous, nothing
        list_publishable = [(choice_PI==1)|((choice_PI==2)&bool_serend) for bool_serend in indices_not_target]

        sources_raw = sources_raw[(sources_raw["EP_EXT_ML"]<6) & indices_not_spurious]
        tab_band_fluxes = [[list(line)] for line in sources_raw["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"]]
        tab_band_fluxerr = []
        nb_src = 0
        nb_alerts = 0
        for line in sources_raw["ERR_EP_1_FLUX","ERR_EP_2_FLUX","ERR_EP_3_FLUX","ERR_EP_4_FLUX","ERR_EP_5_FLUX"]:
            tab_band_fluxerr.append([[list(line)], [list(line)]])
    
        for (ra, dec, pos_err, flux, flux_err, band_flux, band_fluxerr, src_num, angle,
             pn_detml, m1_detml, m2_detml, ep_detml, publishable)\
                in zip(sources_raw["RA"],sources_raw["DEC"],sources_raw["RADEC_ERR"], sources_raw["EP_TOT_FLUX"],\
                        sources_raw["ERR_EP_TOT_FLUX"],tab_band_fluxes, tab_band_fluxerr, sources_raw["SRC_NUM"],
                        off_target_angles, sources_raw["PN_DET_ML"],
                       sources_raw["M1_DET_ML"],sources_raw["M2_DET_ML"],sources_raw["EP_DET_ML"],list_publishable):
            param_holder = ParamHolder()
            param_holder.ra = ra 
            param_holder.dec = dec
            param_holder.pos_err = pos_err
            param_holder.flux = flux
            param_holder.flux_err = flux_err
            param_holder.band_flux = band_flux
            param_holder.band_fluxerr = band_fluxerr
            param_holder.src_num = src_num
            param_holder.off_target_angle = angle
            #param_holder.m1_offax = m1_offax
            #param_holder.m2_offax = m2_offax
            param_holder.pn_detml = pn_detml
            param_holder.m1_detml = m1_detml
            param_holder.m2_detml = m2_detml
            param_holder.ep_detml = ep_detml
            param_holder.choice_pi = choice_PI
            param_holder.publishable = publishable
            param_holder.pointing_type = flag_pointing_type

            print (f"Processing source {src_num}")
            nb_alerts += process_one_source(param_holder, dict_observation_metadata, session, image_data, image_wcs)
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
        
         
def process_one_source(param_holder, observation_metadata, session, image_data, image_wcs):
    tab_alerts=[]
    tab_dic_infos = []
    nb_alerts = 0
    dict_detection_info={}
    dict_detection_info["ObsID"]=observation_metadata["ObsID"]
    dict_detection_info["Date Obs"]=observation_metadata["DateObs"]
    dict_detection_info["Target Name"]=observation_metadata["TargetName"]
    dict_detection_info["Exposure Time"] = str(observation_metadata["ExpTime"])+" s"
    dict_detection_info['SRCNUM'] = str(param_holder.src_num)
    #dict_detection_info['Off-axis Angles'] = f"PN: {param_holder.pn_offax:.1f}', M1: {param_holder.m1_offax:.1f}', M2: {param_holder.m2_offax:.1f}'"
    dict_detection_info['Angular separation from target'] = f"{param_holder.off_target_angle:.1f}'"
    dict_detection_info['Instruments DetML'] = f'PN: {param_holder.pn_detml:.1f}, M1: {param_holder.m1_detml:.1f}, M2: {param_holder.m2_detml:.1f}, EP: {param_holder.ep_detml:.1f}'
    c=SkyCoord(param_holder.ra*u.deg,param_holder.dec*u.deg).to_string('hmsdms', sep=":", precision=1)
    dict_detection_info['Source RA']= f'{np.round(param_holder.ra, 4)}   /   {c.split(" ")[0]}'
    dict_detection_info['Source Dec']= f'{np.round(param_holder.dec, 4)}   /   {c.split(" ")[1]}'
    dict_detection_info['Position Error']=f'{param_holder.pos_err:.2f}"'
    dict_detection_info['ChoicePI']=f'{param_holder.choice_pi}'
    dict_detection_info['Publishable']=f'{param_holder.publishable}'
    dict_detection_info['PointingType']=f'{param_holder.pointing_type}'


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
            ms.save_lightcurve(dict_det_info, flag_alert, image_data, image_wcs)
            ms.save_json_alert(dict_det_info, flag_alert, param_holder)
            nb_alerts += 1
 
    return nb_alerts