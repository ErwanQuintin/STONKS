'''
Created on 28 févr. 2023

@author: michel
'''
import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from tqdm import tqdm
from astropy.constants import c
from LoadSpecificMasterSource import *
import api as rpx
from astropy.coordinates import SkyCoord, search_around_sky, Angle
import time
import shlex
import subprocess

from STONKS_PreComputed_Position_alert import path_to_master_sources, transient_alert

def testing_functions_0804670301():
    """
    Testing function, meant to work on a sub-sample of the catalog corresponding to an observation of NGC 7793 containing
    several variables objects
    """
    #Initial input to get the Ra, Dec, Flux, FluxErr
    pos_err = 1#input("1sigma PosErr of the source (in arcsec)?")
    flux_level = 1e-14#input("Flux Level of the source ?")
    flux_err_level = 0.5e-14#input("Flux Error Level of the source ?")
    band_fluxes = [[1e-12]*5]
    band_fluxerr = [[[3e-13]*5],[[3e-13]*5]]

    src_list_path = os.path.join(path_to_master_sources, "P0804670301EPX000OBSMLI0000.FIT")
    print(f"Loading EPIC source list {src_list_path}")
    raw_data = fits.open(src_list_path, memmap=True)

    #Building Observation information using OBSMLI header
    dict_observation_metadata = {}
    obs_information = raw_data[0].header
    dict_observation_metadata["ObsID"] = str(obs_information['OBS_ID'])
    dict_observation_metadata["Date Obs"] = obs_information['DATE-OBS']
    dict_observation_metadata["Target Name"] = obs_information['OBJECT']
    date = Time(dict_observation_metadata["Date Obs"], format="isot").mjd


    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)

    tab_times = []
    tab_band_fluxes = [[list(line)] for line in sources_raw["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"]]
    tab_band_fluxerr = []
    for line in sources_raw["ERR_EP_1_FLUX","ERR_EP_2_FLUX","ERR_EP_3_FLUX","ERR_EP_4_FLUX","ERR_EP_5_FLUX"]:
        tab_band_fluxerr.append([[list(line)], [list(line)]])

    tab_alerts=[]
    tab_dic_infos = []
    pbar=tqdm(total=len(sources_raw))
    for ra, dec, pos_err, flux, flux_err, band_flux, band_fluxerr, src_num, pn_offax, m1_offax, m2_offax in \
            zip(sources_raw["RA"],sources_raw["DEC"],sources_raw["RADEC_ERR"], sources_raw["EP_TOT_FLUX"],\
                    sources_raw["ERR_EP_TOT_FLUX"],tab_band_fluxes, tab_band_fluxerr, sources_raw["SRC_NUM"],
                sources_raw["PN_OFFAX"],sources_raw["M1_OFFAX"],sources_raw["M2_OFFAX"]):
        start = time.time()
        dict_detection_info={}
        dict_detection_info["ObsID"]=dict_observation_metadata["ObsID"]
        dict_detection_info["Date Obs"]=dict_observation_metadata["Date Obs"]
        dict_detection_info["Target Name"]=dict_observation_metadata["Target Name"]
        dict_detection_info['SRCNUM'] = str(src_num)
        dict_detection_info['Off-axis Angles'] = f'PN: {pn_offax:.1f}", M1: {m1_offax:.1f}", M2: {m2_offax:.1f}"'
        dict_detection_info['Source RA']=np.round(ra, 4)
        dict_detection_info['Source Dec']=np.round(dec, 4)
        dict_detection_info['Position Error']=f'{pos_err:.2f}"'

        result_alert = transient_alert(1, ra, dec, pos_err, flux, flux_err, band_flux,
                                     band_fluxerr, date, var_flag=False)
        if result_alert!=[]: #If there is a transient alert
            tab_alerts += result_alert
            tab_dic_infos.append(dict_detection_info)

        end = time.time()
        tab_times.append(end - start)
        pbar.update(1)
    pbar.close()

    plt.hist(tab_times, bins=np.geomspace(1e-3,1e2,20))
    plt.xscale("log")
    for ms, dict_det_info in zip(tab_alerts,tab_dic_infos):
        ms.save_lightcurve(dict_det_info)
if __name__ == '__main__':
    testing_functions_0804670301()
