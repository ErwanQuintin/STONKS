"""
The Search for Transient Objects in New detections using Known Sources (STONKS) software is meant to allow automatic
comparison between the PPS file of a new XMM-Newton observation, and all the available archival X-ray data. It will send
out an alert and plot the lightcurve and available data for all sources with a long-term variability over a factor of 5.

This requires the script to be run with the same file structure as on the GitHub, to ensure access to hand-made catalog
data. You also need to state the path to the PPS folder.

The STONKS_pipeline_alert.py module uses the archival X-ray catalog loaded through LoadMasterSources.py, compares these
sources to the new detection of a given PPS file, and returns the MasterSource objects for which a variability was
detected.

It is also possible to use the archival catalog strictly for data mining past data, but this is done in the other script
(StudyMasterSources.py).

This software was developed as part of the XMM2ATHENA project. This project has received funding from the European
Union's Horizon 2020 research and innovation programme under grant agreement n°101004168, the XMM2ATHENA project.

Author: Erwan Quintin, erwan.quintin@irap.omp.eu
"""

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
from LoadMasterSources import *
import api as rpx
from astropy.coordinates import SkyCoord, search_around_sky, Angle


#Initial input if you are interested in optical counterparts to X-ray sources, from OM & UVOT data. This is interesting,
#but takes some time to load
input_given = False
while not input_given:
    multiwavelength_information_input = input("Do you want to use OM & UVOT data ? This leads to longer loading time (Please enter Y or N)")
    if multiwavelength_information_input=="Y":
        input_given = True
        multiwavelength_information = True
    elif multiwavelength_information_input=="N":
        input_given = True
        multiwavelength_information = False
    else:
        print("Please answer by Y or N")
print("Choice for OM & UVOT: ", multiwavelength_information)

#This loads all archival data, in the form of a dictionary of Master Sources, each containing several correlated catalog sources
print("\nLoading all archival X-ray data, should take about ~10mn")
dic_master_sources, tab_optical_only = load_master_sources(multiwavelength_information=multiwavelength_information)

#This path leads to two folders, one which corresponds to the PPS data and the other to the Sources data.
path_to_data = input("Absolute path to the folder containing PPS and Sources files?")
#"/home/erwan/Documents/PhD/LongTermVariability/LongTermVariability_Python/STONKS_PipelineVersion/Test_ObsFiles"


def load_new_EP_detections(path_to_data):
    """
    Loads the new EPIC detections from the PPS folder, filters them on quality, and creates the corresponding Source objects
    :parameter
    path_to_data: path to the folder containing PPS data
    :return
    tab_observed_sources: table containing Source objects corresponding to the new detections. See LoadMasterSources.py
    for details on Source and MasterSource objects.
    """
    list_files_pps = os.listdir(f"{path_to_data}/pps/pps/")
    path_to_detection_fits = [f for f in list_files_pps if (("EPX000OBSMLI" in f) and (f.endswith("FTZ")))][0]
    raw_data = fits.open(f"{path_to_data}/pps/pps/"+path_to_detection_fits, memmap=True)
    head = raw_data[0].header
    date = Time(head['DATE-OBS'], format="isot").mjd
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    quality_flags = [elt.count("T")<=3 for elt in sources_raw['EP_FLAG']] #This corresponds to SUM_FLAG<=3
    sources_raw = sources_raw[quality_flags]

    tab_observed_sources = []
    obsid = path_to_detection_fits.split('/')[-1][1:11]
    for XMM_detection in sources_raw:
        current_name = 'NewXMM'+obsid+'-'+str(XMM_detection['SRC_NUM'])
        band_fluxes = [XMM_detection["EP_"+band+"_FLUX"] for band in ("1", "2", "3", "4", "5")]
        band_flux_errors = [XMM_detection["ERR_EP_"+band+"_FLUX"] for band in
                       ("1", "2", "3", "4", "5")]
        current_xmm_source = Source('NewXMM', current_name, [XMM_detection['EP_TOT_FLUX']],
                                    [[XMM_detection['ERR_EP_TOT_FLUX']],[XMM_detection['ERR_EP_TOT_FLUX']]], [date], [band_fluxes], [[band_flux_errors],[band_flux_errors]])
        current_xmm_source.sc_var_flag = False
        current_xmm_source.off_axis = [[XMM_detection["PN_OFFAX"], XMM_detection["M1_OFFAX"], XMM_detection["M2_OFFAX"]]]
        current_xmm_source.var_flags =[False]
        current_xmm_source.ra = XMM_detection['RA']
        current_xmm_source.dec = XMM_detection['DEC']
        current_xmm_source.poserr = 3 * XMM_detection['RADEC_ERR']
        tab_observed_sources.append(current_xmm_source)
    return tab_observed_sources

def load_new_OM_detections():
    """
    Loads the new OM detections from the PPS folder, filters them on quality, and creates the corresponding OpticalSource objects.
    Still WIP, to be completely implemented soon.
    :parameter
    path_to_data: path to the folder containing PPS data
    :return
    tab_observed_sources: table containing OpticalSource objects corresponding to the new detections. See LoadMasterSources.py
    for details on OpticalSource and MasterSource objects.
    """
    list_files_pps = os.listdir(f"{path_to_data}/pps/pps/")
    path_to_OM_detection_fits = [f for f in list_files_pps if (("OMX000OBSMLI" in f) and (f.endswith("FTZ")))][0]
    raw_data = fits.open(path_to_OM_detection_fits, memmap=True)
    head = raw_data[0].header
    date = Time(head['DATE-OBS'], format="isot").mjd
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    #quality_flags = (sources_raw['SUM_FLAG'] <= 2) & (sources_raw['EP_EXTENT'] == 0.0)
    #sources_raw = sources_raw[quality_flags]

    tab_observed_sources = []
    obsid = path_to_detection_fits.split('/')[-1][1:11]
    for XMM_detection in sources_raw:
        current_name = 'NewXMM'+obsid+'_'+str(XMM_detection['SRC_NUM'])
        rates = [XMM_detection["EP_"+band+"_FLUX"]/XMMband_center[band] for band in ("1", "2", "3", "4", "5")]
        rate_errors = [XMM_detection["ERR_EP_"+band+"_FLUX"]/XMMband_center[band] for band in
                       ("1", "2", "3", "4", "5")]
        current_xmm_source = Source('XMM',0, current_name, XMM_detection['EP_TOT_FLUX'],
                                    XMM_detection['ERR_EP_TOT_FLUX'], date, rates, rate_errors)
        current_xmm_source.sc_var_flag = False
        current_xmm_source.off_axis = [[XMM_detection["PN_OFFAX"], XMM_detection["M1_OFFAX"], XMM_detection["M2_OFFAX"]]]
        current_xmm_source.var_flags =[False]
        current_xmm_source.ra = XMM_detection['RA']
        current_xmm_source.dec = XMM_detection['DEC']
        current_xmm_source.poserr = 3 * XMM_detection['RADEC_ERR']
        tab_observed_sources.append(current_xmm_source)
    return tab_observed_sources

def bayesian_likelihood_score(angDist, err1, err2):
    """
    Likelihood score of a given match, computed by Bayes formula and convoluting two 2D gaussians. There is no prior in
    this formula, opposed to NWAY, because this will be used to compare two scores between matches in similar catalogs,
    so priors cancel out.
    :param angDist: Angular distance between sky positions of the two sources
    :param err1: 2D 1sigma position error of the first source
    :param err2: 2D 1sigma position error of the second source
    :return: the likelihood score of this particular match, to be compared to other competing matches
    """
    return (2 / (err1**2 + err2**2)) * np.exp(-angDist / (err1**2 + err2**2))

def match_pps_with_MasterSources(tab_new_sources):
    """
    Matches the position of the new sources with the archival data, keeping the possible ambiguous matches up to 20".
    :param tab_new_sources: table containing the new detected Source objects, output from load_new_EP_detections()
    :return: dic_splitted_associations: Dictionnary containing for keys an identifier of each new source, and as values
    a table containing all possible matches for each source. In this table, each match is itself a table, containing the
    identifiers of both the new Source and the corresponding MasterSource, their separation, and the likelihood score of
    the match computed using bayesian_likelihood_score()
    """
    tab_new_ra = [src.ra for src in tab_new_sources]
    tab_new_dec = [src.dec for src in tab_new_sources]
    tab_new_poserr = [src.poserr for src in tab_new_sources]
    tab_new_srcnum = [src.name.split("-")[-1] for src in tab_new_sources]
    new_positions = SkyCoord(tab_new_ra, tab_new_dec, unit='deg')
    tab_archival_ra, tab_archival_dec, tab_archival_poserr, tab_archival_ms_id = [], [], [], []
    for ms in tqdm(dic_master_sources.values()):
        tab_archival_ra.append(ms.ra)
        tab_archival_dec.append(ms.dec)
        tab_archival_poserr.append(3*ms.pos_err)
        tab_archival_ms_id.append(ms.id)
    archival_positions = SkyCoord(tab_archival_ra, tab_archival_dec, unit='deg')

    #At this we have positions of the new sources and the archival sources. We now match them up to 20",
    #regardless of the actual position errors
    max_err = Angle(20, unit="arcsec")
    idx1, idx2, sep2D, sep3D = search_around_sky(new_positions, archival_positions, max_err)
    good_ind1 = []
    good_ind2 = []
    good_sep = []
    for ind1, ind2, sep in zip(idx1, idx2, sep2D):
        #We now only keep matches for which the separation is below the sum of 3 sigma error bars
        if sep.arcsec < max(tab_new_poserr[ind1], 1) + max(tab_archival_poserr[ind2], 1):
            good_ind1.append(ind1)
            good_ind2.append(ind2)
            good_sep.append(sep)

    #We then split the table containing all associations in sub-tables, each containing all associations for one new Source
    dic_splitted_associations = {}
    for ind1, ind2, sep in zip(good_ind1, good_ind2, good_sep):
        if ind1 in dic_splitted_associations.keys():
            dic_splitted_associations[ind1].append([ind1, tab_archival_ms_id[ind2], sep, bayesian_likelihood_score(sep.arcsec, tab_new_poserr[ind1]/3,tab_archival_poserr[ind2]/3)])
        else:
            dic_splitted_associations[ind1] = [[ind1, tab_archival_ms_id[ind2], sep, bayesian_likelihood_score(sep.arcsec, tab_new_poserr[ind1]/3,tab_archival_poserr[ind2]/3)]]
    return dic_splitted_associations

def solve_associations_new_detections(tab_new_sources, dic_splitted_associations):
    """
    Takes the dictionary of possible associations and decides on the ambiguous cases
    :param tab_new_sources: Table containing the Source objects for new detections
    :param dic_splitted_associations: Dictionary containing all possible matches for each of these new Sources, if any
    :return: tab_new_sources_for_new_ms: Sub-table of tab_new_sources with only Source objects with no archival match
    tab_new_sources_for_archival_ms: Sub_table of tab_new_sources with only Source objects matching unambiguously archival data
    tab_indices_impacted_master_sources: Indices of the archival MasterSource objects corresponding to these matches
    """
    tab_indices_impacted_master_sources = []
    tab_new_sources_for_archival_ms = []
    tab_new_sources_for_new_ms = []

    count_archival_unambiguous = 0
    count_archival_unambiguousXMM = 0
    count_archival_unambiguousBayesian = 0
    count_archival_ambiguous = 0
    count_new_sources = 0
    for new_source_ind, new_source in enumerate(tab_new_sources):
        if new_source_ind in dic_splitted_associations.keys():
            #It matched with the catalog
            if len(dic_splitted_associations[new_source_ind])>1:
                #It was an ambiguous matching: start by checking if any archival XMM, then using Bayesian match comparisons
                number_of_archival_XMM_matches = len([association[1] for association in dic_splitted_associations[new_source_ind] if "XMM" in dic_master_sources[association[1]].sources.keys()])
                if number_of_archival_XMM_matches == 1:
                    #Single XMM match takes priority
                    tab_indices_impacted_master_sources.append([association[1] for association in dic_splitted_associations[new_source_ind] if "XMM" in dic_master_sources[association[1]].sources.keys()][0])
                    tab_new_sources_for_archival_ms.append(new_source)
                    count_archival_unambiguousXMM+=1
                else:
                    #We need Bayesian decision framework to treat this association: over a factor 5 difference in the scores between competing matches
                    all_B_factors = [association[-1] for association in dic_splitted_associations[new_source_ind]]
                    ordered_B = np.argsort(all_B_factors)[::-1]
                    if all_B_factors[ordered_B[0]]>5*all_B_factors[ordered_B[1]]:
                        #Not that ambiguous
                        tab_indices_impacted_master_sources.append(dic_splitted_associations[new_source_ind][ordered_B[0]][1])
                        tab_new_sources_for_archival_ms.append(new_source)
                        count_archival_unambiguousBayesian+=1
                    else:
                        #Too ambiguous to decide
                        count_archival_ambiguous+=1
            else:
                #It was an unambiguous matching
                count_archival_unambiguous += 1
                tab_indices_impacted_master_sources.append(dic_splitted_associations[new_source_ind][0][1])
                tab_new_sources_for_archival_ms.append(new_source)
        else:
            #It's a new source
            tab_new_sources_for_new_ms.append(new_source)
            count_new_sources+=1
    print(count_new_sources, count_archival_unambiguous, count_archival_unambiguousXMM, count_archival_unambiguousBayesian, count_archival_ambiguous)
    return tab_new_sources_for_new_ms, tab_new_sources_for_archival_ms, tab_indices_impacted_master_sources

def compute_upper_limit(ms):
    """
    Computes the XMM-Newton upper limit on the position of a given MasterSource object, using the RapidXMM framework from
    Ruiz et al. 21, available in api.py
    At this date, there is an issue with ESAC servers so this does not work. Should be fixed soon.
    :param ms: the MasterSource object, corresponding to a new detection, on which to compute the Upper Limits
    :return: xmm_ul: table containing the 3 sigma band 8 EPIC flux upper limits
        xmm_ul_dates: table containing the dates of the upper limits
        xmm_ul_obsids: table containing the ObsIDs of the upper limits
    """
    #Dictionary containing the count-rates to flux conversion factors, in various instruments and filter, in band 8
    ecf_dic = {"PN":{"Open":4.1236,"Thin1":3.3243,"Thin2":3.3243,"Medium":3.1924,"Thick":2.5928},
               "M1":{"Open":1.0916,"Thin1":0.929,"Thin2":0.929,"Medium":0.900,"Thick":0.7779},
               "M2":{"Open":1.1003,"Thin1":0.9358,"Thin2":0.9358,"Medium":0.906,"Thick":0.7829}}
    #tab_upper_limits = rpx.query_coords(SkyCoord(ra=ms.ra, dec=ms.dec, unit="deg"))
    tab_upper_limits = rpx.query_coords([SkyCoord(ra=184.45, dec=+29.9, unit="deg")])
    print(tab_upper_limits)
    #tab_ecf = [ecf_dic[line["instrum"]][line["filt"]]*1e11 for line in sources_raw]
    #tab_ul_flux8_3sig = [ul_rate/ecf for ecf, ul_rate in zip(tab_ecf, sources_raw['band8_ul_sigma3'])]

def treat_new_detections(tab_new_sources_for_new_ms, tab_new_sources_for_archival_ms, tab_indices_impacted_master_sources):
    """
    Takes the various Source objects divided between those known in archives and the new ones. For the former, adds the
    new detection treated as an additional Source object. For the latter, creates a new MasterSource object for it and
    computes the archival XMM-Newton upper limits, if any.
    :param tab_new_sources_for_new_ms: Sub-table of tab_new_sources with only Source objects with no archival match
    :param tab_new_sources_for_archival_ms: Sub_table of tab_new_sources with only Source objects matching unambiguously archival data
    :param tab_indices_impacted_master_sources: Identifiers of the archival MasterSource objects corresponding to these matches
    :return: master_sources_for_transient_alert: a Table containing all the MasterSource objects concerned by this new observation
    """
    last_ms_id = max(list(dic_master_sources.keys()))
    master_sources_for_transient_alert = []
    for ind_new_source, new_source in enumerate(tab_new_sources_for_new_ms):
        ms = MasterSource(last_ms_id+ind_new_source+1, [new_source], new_source.ra, new_source.dec, new_source.poserr, [])
        #ul = compute_ul(ms)
        master_sources_for_transient_alert.append(ms)

    for ind_ms, new_source in zip(tab_indices_impacted_master_sources, tab_new_sources_for_archival_ms):
        ms = dic_master_sources[ind_ms]
        new_ms = MasterSource(ms.id,  list(ms.sources.values())+[new_source], ms.ra, ms.dec, ms.pos_err, list(ms.optical_sources.values()))
        dic_master_sources[ind_ms] = new_ms
        master_sources_for_transient_alert.append(new_ms)
    return master_sources_for_transient_alert

def transient_alert(path_to_data):
    """
    Sends out alerts in the case of a transient object
    :param path_to_data: Path to the folder containing the PPS files
    :return: tab_alerts: Sub-table of the MasterSource objects variable enough to trigger an alert. We also plot their lightcurves
    """
    #Load the new detections
    tab_new_sources = load_new_EP_detections(path_to_data)

    #Matches them with archival data
    dic_splitted_associations = match_pps_with_MasterSources(tab_new_sources)

    #Solves the ambiguous matches
    tab_new_sources_for_new_ms, tab_new_sources_for_archival_ms, tab_indices_impacted_master_sources = solve_associations_new_detections(
        tab_new_sources, dic_splitted_associations)

    #Sums up the MasterSource objects concerned by this new detection
    master_sources_for_transient_alert = treat_new_detections(tab_new_sources_for_new_ms,
                                                              tab_new_sources_for_archival_ms,
                                                              tab_indices_impacted_master_sources)

    #Check if any of these MasterSource objects is variable
    tab_alerts = []
    for ms in master_sources_for_transient_alert:
        if ms.var_ratio > 5:
            tab_alerts.append(ms)
    print("Number of alerts:", len(tab_alerts))
    for ms in tab_alerts:
        ms.plot_lightcurve()
    return tab_alerts


tab_alerts = transient_alert(path_to_data)