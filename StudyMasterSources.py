import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from tqdm import tqdm
from mocpy import MOC
from astropy.coordinates import SkyCoord
import webbrowser
import arviz as az
import seaborn as sns
from scipy.stats import pearsonr
from matplotlib.colors import ListedColormap, LogNorm
from itertools import combinations
from scipy.stats import gaussian_kde, binned_statistic_2d, linregress, binned_statistic
from astropy.constants import c
from matplotlib.gridspec import GridSpec
import imageio
from PIL import Image
import os
colormap = ListedColormap(np.loadtxt('colormap.txt')/255.)
import cmasher as cmr

from LoadMasterSources import load_master_sources, catalogs, optical_catalogs

dic_master_sources, tab_optical_only = load_master_sources(multiwavelength_information=True)
optical_bands = ["UVW2", "UVM2", "UVW1", "U", "B", "V"]


def number_of_alerts(start, end):
    #colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    count_over_threshold3 = []
    count_over_threshold3_errorbar = []
    count_over_threshold5 = []
    count_over_threshold5_errorbar = []
    dates = []
    nbr_steps = int((end-51577)/60)
    for step in tqdm(range(nbr_steps)):
        date_limit_low = start + (end-start)*(step/nbr_steps)
        dates.append(date_limit_low)
        date_limit_up = date_limit_low + 60
        dic_archival_fluxes_upper, dic_archival_fluxes_lower, dic_new_fluxes_upper, dic_new_fluxes_lower, dic_new_times = separate_detections_dates(date_limit_low, date_limit_up)
        tab_alert_times_3 = []
        tab_alert_times_5 = []
        for ms_id in dic_master_sources.keys():
            if len(dic_archival_fluxes_upper[ms_id])>0 and len(dic_new_fluxes_upper[ms_id])>0:
                if (np.nanmin(dic_new_fluxes_upper[ms_id])<np.nanmin(dic_archival_fluxes_upper[ms_id])):
                    if (np.nanmin(dic_new_fluxes_upper[ms_id])<(1/3)*np.nanmax(dic_archival_fluxes_lower[ms_id])):
                        time=dic_new_times[ms_id][np.argmin(dic_new_fluxes_upper[ms_id])]
                        tab_alert_times_3.append(time)
                        if (np.nanmin(dic_new_fluxes_upper[ms_id])<(1/5)*np.nanmax(dic_archival_fluxes_lower[ms_id])):
                            tab_alert_times_5.append(time)
                elif (np.nanmax(dic_new_fluxes_lower[ms_id])>np.nanmax(dic_archival_fluxes_lower[ms_id])):
                    if (np.nanmax(dic_new_fluxes_lower[ms_id])>5*np.nanmin(dic_archival_fluxes_upper[ms_id])):
                        time = dic_new_times[ms_id][np.argmax(dic_new_fluxes_lower[ms_id])]
                        tab_alert_times_3.append(time)
                        if (np.nanmax(dic_new_fluxes_lower[ms_id])>5*np.nanmin(dic_archival_fluxes_upper[ms_id])):
                            tab_alert_times_5.append(time)

        count_over_threshold3.append(len(tab_alert_times_3)/60)
        count_over_threshold3_errorbar.append(np.std(np.histogram(tab_alert_times_3, bins=60)[0]))
        count_over_threshold5.append(len(tab_alert_times_5) / 60)
        count_over_threshold5_errorbar.append(np.std(np.histogram(tab_alert_times_5, bins=60)[0]))
    #plt.plot(dates, count_over_threshold3, label="Var > 3")
    #plt.fill_between(dates, np.array(count_over_threshold3)-np.array(count_over_threshold3_errorbar), np.array(count_over_threshold3)+np.array(count_over_threshold3_errorbar),alpha=0.25)
    #plt.axhline(y=np.nanmean(count_over_threshold3), ls="--",c=colors[0])
    plt.plot(dates, count_over_threshold5, label="Var > 5")
    plt.fill_between(dates, np.array(count_over_threshold5) - np.array(count_over_threshold5_errorbar),
                     np.array(count_over_threshold5) + np.array(count_over_threshold5_errorbar), alpha=0.25)
    plt.axhline(y=np.nanmean(count_over_threshold5), ls="--")#,c=colors[1])
    plt.xlabel("Starting date of 2 months alert test")
    plt.ylabel("Number of alerts per day (2-months average)")
    plt.legend()
    plt.show()
number_of_alerts(51577, 59200)