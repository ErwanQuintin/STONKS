"""
The Search for Transient Objects in New detections using Known Sources (STONKS) software is meant to allow automatic
comparison between the PPS file of a new XMM-Newton observation, and all the available archival X-ray data. It will send
out an alert and plot the lightcurve and available data for all sources with a long-term variability over a factor of 5.

The STONKS_PrecComputed_Position_alert.py module uses the archival X-ray catalog loaded through LoadSpecificMasterSource.py,
 compares the sources around the target position to the new detection at a given flux level, and returns the
 long-term lightcurve of the associated MasterSource object and an alert if it is variable.

It is also possible to use the archival catalog strictly for data mining past data, but this is done in the other script
(StudyMasterSources.py).

This software was developed as part of the XMM2ATHENA project. This project has received funding from the European
Union's Horizon 2020 research and innovation programme under grant agreement n°101004168, the XMM2ATHENA project.

Author: Erwan Quintin, erwan.quintin@irap.omp.eu
"""

import os
import numpy as np
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from tqdm import tqdm
from astropy.constants import c
from core.LoadSpecificMasterSource import MasterSource, Source, load_specific_master_sources
from core import api as rpx
from astropy.coordinates import SkyCoord, search_around_sky, Angle
import time
import shlex
import subprocess
from constants import PATHTO
from scipy.special import gammaincinv, gammaincc, gammainc
from astroquery.simbad import Simbad
Simbad.add_votable_fields("otype")



def load_master_sources_positions(obsid, ra_target, dec_target):
    """
    This will build the sub-table of the MasterSource catalog corresponding to 30' radius zone, if it doesn't already exist.
    It then loads this archival data, to be then compared to the new detection.
    :param obsid: ID of the Observation, used to check if the sub-region of the catalog was already computed
    :param ra_target: RA of the source
    :param dec_target: Dec of the source
    :return: Identification, positional and flux information of the MasterSources in the 45' radius.
    """
    list_precomputed_obsids = os.listdir(os.path.join(PATHTO.master_sources,'PreComputedObsidMatches'))
    list_precomputed_obsids=[elt.split(".")[0] for elt in list_precomputed_obsids]
    if str(obsid) not in list_precomputed_obsids or f"UpperLimits_{obsid}" not in list_precomputed_obsids:
        cmd = f"{PATHTO.stilts_cmd} tpipe {os.path.join(PATHTO.master_sources,'Master_source_HistoricalExtremes.fits')} cmd='select \"skyDistanceDegrees(MS_RA,MS_DEC,{ra_target},{dec_target})*60<45\"' \
        out={os.path.join(PATHTO.precomputed_obsids, str(obsid)+'.fits')}"
        cmd = shlex.split(cmd)
        subprocess.run(cmd)

        cmd = f"{PATHTO.stilts_cmd} tpipe {os.path.join(PATHTO.master_sources,'Master_source_XMM_UpperLimits.fits')} cmd='select \"skyDistanceDegrees(MS_RA,MS_DEC,{ra_target},{dec_target})*60<45\"' \
                out={os.path.join(PATHTO.precomputed_obsids, 'UpperLimits_'+str(obsid)+'.fits')}"
        cmd = shlex.split(cmd)
        subprocess.run(cmd)

    raw_data = fits.open(f"{os.path.join(PATHTO.master_sources,'PreComputedObsidMatches',str(obsid)+'.fits')}", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    tab_positions = SkyCoord(sources_raw["MS_RA"],sources_raw["MS_DEC"], unit='deg')
    tab_poserr = 3*sources_raw["MS_POSERR"]
    tab_has_xmm = sources_raw["XMM"]!=' '
    return sources_raw["MS_ID"], tab_positions, tab_poserr,tab_has_xmm, sources_raw['HistoricalMin'],sources_raw['HistoricalMax'],sources_raw['LastSeenFlux'],sources_raw['LastSeenFluxErr'],sources_raw['LowerUpperLimit']

def create_new_Source(ra_target, dec_target, pos_err, flux, flux_err, band_fluxes, band_fluxerr, date):
    """
    Takes all the information of the new detection and creates a Source object with it
    """
    name="New Detection"
    new_source = Source('NewXMM', name, [flux],
           [[flux_err], [flux_err]], [date], band_fluxes,  band_fluxerr)
    new_source.ra = ra_target
    new_source.dec = dec_target
    new_source.poserr = 3 * pos_err
    return new_source

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

def match_newDet_with_MasterSources_position(tab_ms_pos, tab_ms_poserr, new_source_pos, new_source_poserr):
    """
    Matches the position of the new sources with the archival data, keeping the possible ambiguous matches up to 20".
    :param tab_new_sources: table containing the new detected Source objects, output from load_new_EP_detections()
    :return: dic_splitted_associations: Dictionnary containing for keys an identifier of each new source, and as values
    a table containing all possible matches for each source. In this table, each match is itself a table, containing the
    identifiers of both the new Source and the corresponding MasterSource, their separation, and the likelihood score of
    the match computed using bayesian_likelihood_score()
    """

    #At this point we have positions of the new sources and the archival sources. We now match them up to 20",
    #regardless of the actual position errors
    max_err = Angle(20, unit="arcsec")
    idx1, idx2, sep2D, sep3D = search_around_sky(new_source_pos, tab_ms_pos, max_err)
    good_ind1 = []
    good_ind2 = []
    good_sep = []
    for ind1, ind2, sep in zip(idx1, idx2, sep2D):
        #We now only keep matches for which the separation is below the sum of 3 sigma error bars
        if sep.arcsec < max(3*new_source_poserr, 1) + max(tab_ms_poserr[ind2], 1):
            good_ind1.append(ind1)
            good_ind2.append(ind2)
            good_sep.append(sep)

    #We then split the table containing all associations in sub-tables, each containing all associations for one new Source
    list_possible_associations = []
    for ind1, ind2, sep in zip(good_ind1, good_ind2, good_sep):
        list_possible_associations.append([ind1, ind2, sep, bayesian_likelihood_score(sep.arcsec, new_source_poserr/3,tab_ms_poserr[ind2]/3)])
    return list_possible_associations

def solve_associations_new_detections(tab_has_xmm, list_possible_associations):
    """
    Takes the dictionary of possible associations and decides on the ambiguous cases
    :param tab_new_sources: Table containing the Source objects for new detections
    :param dic_splitted_associations: Dictionary containing all possible matches for each of these new Sources, if any
    :return: tab_new_sources_for_new_ms: Sub-table of tab_new_sources with only Source objects with no archival match
    tab_new_sources_for_archival_ms: Sub_table of tab_new_sources with only Source objects matching unambiguously archival data
    tab_indices_impacted_master_sources: Indices of the archival MasterSource objects corresponding to these matches
    """
    index_impacted_master_sources = []
    if len(list_possible_associations)!=0:
        #It matched with the catalog
        if len(list_possible_associations)>1:
            #It was an ambiguous matching: start by checking if any archival XMM, then using Bayesian match comparisons
            number_of_archival_XMM_matches = np.sum([tab_has_xmm[association[1]] for association in list_possible_associations])
            if number_of_archival_XMM_matches == 1:
                #Single XMM match takes priority
                index_impacted_master_sources.append([association[1] for association in list_possible_associations if tab_has_xmm[association[1]]])
                flag_known_ms=1
            else:
                #We need Bayesian decision framework to treat this association: over a factor 5 difference in the scores between competing matches
                all_B_factors = [association[-1] for association in list_possible_associations]
                ordered_B = np.argsort(all_B_factors)[::-1]
                if all_B_factors[ordered_B[0]]>5*all_B_factors[ordered_B[1]]:
                    #Not that ambiguous
                    index_impacted_master_sources.append(list_possible_associations[ordered_B[0]][1])
                    flag_known_ms=1
                else:
                    flag_known_ms = -1

        else:
            #It was an unambiguous matching
            flag_known_ms = 1
            index_impacted_master_sources.append(list_possible_associations[0][1])
    else:
        #It's a new source
        flag_known_ms = 0
    return flag_known_ms, index_impacted_master_sources

def ul_from_counts_bck(N, B, cl):
    return gammaincinv(N+1, gammainc(N+1,B)+cl*gammaincc(N+1,B))-B

def combine_EPIC_ul(tab_pointed_upper_limits, ecf_dic):
    """Takes the ordered table of count upper limits from all EPIC instruments for each observations, combines them and
    converts them to flux"""
    tab_splitted_ul_info = np.split(tab_pointed_upper_limits, np.unique(tab_pointed_upper_limits["obsid"], return_index=True)[1][1:])
    xmm_ul=[]
    xmm_ul_dates=[]
    for observation in tab_splitted_ul_info:
        N_EPIC = np.sum([line[24] for line in observation])
        B_EPIC = np.sum([line[25] for line in observation])
        lower_fraction = np.sum([line[23]*line[9]*ecf_dic[line[30]][line[31].strip()]*1e11 for line in observation])
        xmm_ul.append(ul_from_counts_bck(N_EPIC, B_EPIC, 0.997)/lower_fraction)
        xmm_ul_dates.append(observation[0][7])
    xmm_ul_dates=Time(xmm_ul_dates,format="isot").mjd
    return xmm_ul, xmm_ul_dates

def compute_upper_limit(ra, dec, flux, date_obs):
    """
    Computes the XMM-Newton upper limit on a given position, using the RapidXMM framework from
    Ruiz et al. 21, available in api.py
    :param ra, dec: the position of the object, corresponding to a new detection, on which to compute the Upper Limits
    :return: xmm_ul: table containing the 3 sigma band 8 EPIC flux upper limits
        xmm_ul_dates: table containing the dates of the upper limits
        xmm_ul_obsids: table containing the ObsIDs of the upper limits
    """
    #Dictionary containing the count-rates to flux conversion factors, in various instruments and filter, in band 8
    ecf_dic = {"PN":{"Open":4.1236,"Thin1":3.3243,"Thin2":3.3243,"Medium":3.1924,"Thick":2.5928},
               "M1":{"Open":1.0916,"Thin1":0.929,"Thin2":0.929,"Medium":0.900,"Thick":0.7779},
               "M2":{"Open":1.1003,"Thin1":0.9358,"Thin2":0.9358,"Medium":0.906,"Thick":0.7829}}

    if flux>1e-12:
        tab_upper_limits = rpx.query_radec(ra=[ra], dec=[dec])
    else:
        tab_upper_limits = rpx.query_radec(ra=[ra], dec=[dec], obstype="pointed")
    xmm_ul, xmm_ul_dates, slew_ul, slew_ul_dates= [],[],[],[]
    if len(tab_upper_limits)>0:
        ul_dates = Time(tab_upper_limits['start_date'],format="isot").mjd
        tab_pointed_upper_limits = tab_upper_limits[np.where((tab_upper_limits["obstype"]=="pointed") & (tab_upper_limits["band8_flags"]==0) & (ul_dates<date))]
        if len(tab_pointed_upper_limits)>0:
            tab_pointed_upper_limits=tab_pointed_upper_limits[np.argsort(tab_pointed_upper_limits["obsid"])]
            xmm_ul,xmm_ul_dates = combine_EPIC_ul(tab_pointed_upper_limits, ecf_dic)#[ul_rate/ecf for ecf, ul_rate in zip(tab_ecf, tab_pointed_upper_limits['band8_ul_sigma3'])]

        tab_slew_upper_limits = tab_upper_limits[np.where((tab_upper_limits["obstype"]=="slew") & (tab_upper_limits["band8_flags"]==0) & (ul_dates<date))]
        tab_ecf = [ecf_dic[line["instrum"]][line["filt"]] * 1e11 for line in tab_slew_upper_limits]
        slew_ul = [ul_rate / ecf for ecf, ul_rate in zip(tab_ecf, tab_slew_upper_limits['band8_ul_sigma3'])]
        slew_ul_dates = Time(tab_slew_upper_limits['start_date'],format="isot").mjd
    return xmm_ul, xmm_ul_dates, slew_ul, slew_ul_dates

def match_Simbad(ra_target, dec_target, pos_err):
    dic_classifier = {'': '',
                      'X': 'Unknown', 'IR': 'Unknown', 'Radio': 'Unknown', 'MIR': 'Unknown',
                      'NIR': 'Unknown', 'HH': 'Unknown', 'HI': 'Unknown', 'HII': 'Unknown',
                      'LensedImage': 'Unknown', 'LensingEv': 'Unknown', 'Maser': 'Unknown',
                      'MolCld': 'Unknown', 'PartofCloud': 'Unknown', 'Radio(sub-mm)': 'Unknown',
                      'Blue': 'Unknown', 'Possible_lensImage': 'Unknown', 'Unknown': 'Unknown',
                      'Radio(mm)': 'Unknown', 'denseCore': 'Unknown', 'Radio(cm)': 'Unknown',
                      'UV': 'Unknown', 'PN': 'Unknown', 'PN?': 'Unknown', "EmObj": 'Unknown',
                      'DkNeb': 'Unknown', 'Transient': 'Unknown', 'Candidate_LensSystem': 'Unknown',
                      'FIR': 'Unknown', 'multiple_object': 'Unknown', 'GravLensSystem': 'Unknown',
                      'Bubble': 'Unknown', 'Cloud': 'Unknown', 'SFregion': 'Unknown',
                      'Inexistent': 'Unknown', 'gamma': 'Unknown', 'GravLens': 'Unknown',
                      'HVCld': 'Unknown', 'Candidate_Lens': 'Unknown', 'ISM': 'Unknown',
                      'Void': 'Unknown', 'RfNeb': 'Unknown', 'HIshell': 'Unknown', 'Outflow': 'Unknown',
                      'radioBurst': 'Unknown', 'Region': 'Unknown', 'Globule': 'Unknown',
                      'outflow?': 'Unknown', 'ComGlob': 'Unknown',
                      'GinCl': 'Galaxy', 'Galaxy': 'Galaxy', 'AGN': 'Galaxy', 'GiC': 'Galaxy', 'Sy1': 'Galaxy',
                      'Sy2': 'Galaxy', 'AGN_Candidate': 'Galaxy', 'QSO': 'Galaxy', 'Seyfert_1': 'Galaxy',
                      'Seyfert_2': 'Galaxy', 'LINER': 'Galaxy', 'EmG': 'Galaxy', 'RadioG': 'Galaxy', 'BClG': 'Galaxy',
                      'LSB_G': 'Galaxy', 'LensedG': 'Galaxy', 'LensedQ': 'Galaxy', 'GroupG': 'Galaxy',
                      'PartOfG': 'Galaxy', 'BLLac': 'Galaxy', 'GinPair': 'Galaxy', 'Possible_ClG': 'Galaxy',
                      'Possible_G': 'Galaxy', 'Possible_GrG': 'Galaxy', 'GinGroup': 'Galaxy', 'HII_G': 'Galaxy',
                      'Blazar': 'Galaxy', 'ClG': 'Galaxy', 'QSO_Candidate': 'Galaxy', 'Seyfert': 'Galaxy',
                      'Blazar_Candidate': 'Galaxy', 'StarburstG': 'Galaxy', 'IG': 'Galaxy', 'SuperClG': 'Galaxy',
                      'PartofG': 'Galaxy', 'Compact_Gr_G': 'Galaxy', 'PairG': 'Galaxy', 'BLLac_Candidate': 'Galaxy',
                      'BlueCompG': 'Galaxy',
                      'Orion_V*': 'Star', 'TTau*': 'Star', 'EB*': 'Star', 'YSO': 'Star', 'SB*': 'Star', '**': 'Star',
                      'Star': 'Star', 'RotV*': 'Star', 'Candidate_RGB*': 'Star', 'low-mass*': 'Star', 'V*': 'Star',
                      'PulsV*': 'Star', 'AGB*': 'Star', 'S*': 'Star', 'Candidate_YSO': 'Star', 'PM*': 'Star',
                      'Irregular_V*': 'Star', 'Em*': 'Star', 'LPV*': 'Star', 'Mira': 'Star', 'WR*': 'Star',
                      'Pec*': 'Star', 'Planet?': 'Star', 'Planet': 'Star', 'Eruptive*': 'Star', 'Cl*': 'Star',
                      'OpCl': 'Star', 'Assoc*': 'Star', 'PulsV*WVir': 'Star', 'PulsV*bCep': 'Star', 'RRLyr': 'Star',
                      'C*': 'Star', 'EllipVar': 'Star', 'Candidate_EB*': 'Star', 'Candidate_PulsV*WVir': 'Star',
                      'Candidate_LP*': 'Star', 'pulsV*SX': 'Star', 'Candidate_RSG*': 'Star', 'BYDra': 'Star',
                      'Be*': 'Star', 'Candidate_RRLyr': 'Star', 'BlueSG*': 'Star', 'Erupt*RCrB': 'Star', 'RGB*': 'Star',
                      'RSCVn': 'Star', 'gammaDor': 'Star', 'Cl*?': 'Star', 'Candidate_C*': 'Star', 'HB*': 'Star',
                      'Cepheid': 'Star', 'Ae*': 'Star', 'Candidate_TTau*': 'Star', 'deltaCep': 'Star',
                      'HotSubdwarf': 'Star', 'Candidate_AGB*': 'Star', 'YellowSG*': 'Star', 'Symbiotic*': 'Star',
                      'PulsV*delSct': 'Star', 'BlueStraggler': 'Star', 'Candidate_post-AGB*': 'Star',
                      'RotV*alf2CVn': 'Star', 'OH/IR': 'Star', 'V*?': 'Star', 'Candidate_BSG*': 'Star',
                      'RedSG*': 'Star', 'Candidate_brownD*': 'Star', 'Candidate_Mi*': 'Star', 'Candidate_HB*': 'Star',
                      'Candidate_Be*': 'Star', 'Candidate_SN*': 'Star', 'brownD*': 'Star', 'SG*': 'Star',
                      'PulsV*RVTau': 'Star', 'Candidate_WR*': 'Star', 'HV*': 'Star', 'Candidate_Hsd': 'Star',
                      'Candidate_Ae*': 'Star', 'Candidate_Cepheid': 'Star', 'post-AGB*': 'Star', 'Candidate_**': 'Star',
                      'Candidate_Symb*': 'Star', 'Candidate_S*': 'Star', 'Candidate_SG*': 'Star',
                      'Candidate_low-mass*': 'Star',
                      'GlCl': 'CompactObject', 'GlCl?': 'CompactObject', 'Pulsar': 'CompactObject',
                      'ULX': 'CompactObject', 'ULX?': 'CompactObject', 'HMXB': 'CompactObject',
                      'Candidate_HMXB': 'CompactObject', 'LMXB': 'CompactObject', 'Candidate_LMXB': 'CompactObject',
                      'Nova': 'CompactObject', 'CataclyV*': 'CompactObject', 'XB': 'CompactObject',
                      'SNR': 'CompactObject', 'SNR?': 'CompactObject', 'Candidate_WD*': 'CompactObject',
                      'WD*': 'CompactObject', 'SN': 'CompactObject', 'gammaBurst': 'CompactObject',
                      'Candidate_XB*': 'CompactObject', 'Candidate_BH': 'CompactObject', 'NS': 'CompactObject',
                      'Candidate_NS': 'CompactObject', 'Neutron*': 'CompactObject', 'Candidate_CV*': 'CompactObject',
                      'Candidate_Nova': 'CompactObject','HighMassXBin':'CompactObject'}

    result_table = Simbad.query_region(SkyCoord(ra_target, dec_target,
                                                      unit=(u.deg, u.deg), frame='icrs'),
                                       radius=10*u.arcsec)
    # Id result_table is none Simbad raises a script error (don't mind)
    if not result_table is None:
        result = result_table[0]
        print(result.keys())
        if result["otype"] in dic_classifier.keys():
            simbad_type = dic_classifier[result["otype"]]
        else:
            simbad_type = "Unknown"
        simbad_name = result["main_id"]
    else:
        simbad_type = ""
        simbad_name = ""
    return simbad_type, simbad_name

def transient_alert(session, obsid, ra_target, dec_target, pos_err, flux, flux_err, band_fluxes, band_fluxerr, date, var_flag, ul=True):
    """
    Sends out alerts in the case of a transient object
    :param obsid, ra_target, dec_target, pos_err, flux, flux_err, band_fluxes, band_fluxerr, date, src_num, var_flag:
    all information relative to the new detection. Band_fluxes corresponds to the fluxes in the 5 bands.
    src_num is the identifier from the OBSMLI file, used as filename for the saved lightcurve.
    var_flag is the flag for short-term variability.
    ul is a boolean, if True then the upper limits are computed in the case of a new MasterSource.
    :return: tab_alerts: list containing the MasterSource if it's variable (in which case we save a PDF of the
    lightcurve as well), or empty list if not variable.
    """
    # This loads all relevant archival data, in the form of a dictionary of Master Sources, each containing several correlated catalog sources
    tab_ms_id, tab_ms_pos, tab_ms_poserr,tab_has_xmm,tab_min, tab_max, tab_last, tab_last_err, tab_ul = load_master_sources_positions(obsid,ra_target, dec_target)

    #Matches them with archival data
    new_source_pos = SkyCoord([ra_target],[dec_target],unit="deg")
    list_possible_associations = match_newDet_with_MasterSources_position(tab_ms_pos, tab_ms_poserr, new_source_pos,pos_err)

    #Solves the ambiguous matches
    flag_known_ms, index_impacted_master_sources = solve_associations_new_detections(tab_has_xmm, list_possible_associations)

    tab_alerts=[]
    flag_alerts=[]
    info_source=[]
    if flag_known_ms==1:
        #flux=tab_last[0] #Temporary test, should send alerts only if already archival variable
        mini_hist = np.nanmin((tab_min[index_impacted_master_sources[0]],tab_ul[index_impacted_master_sources[0]]))
        if flux < flux_err:
            low_flux = 0
        else:
            low_flux = flux-flux_err
        var_ratio = np.nanmax((low_flux/mini_hist, tab_max[index_impacted_master_sources[0]]/(flux+flux_err)))
        if (var_ratio > 5) :#or var_flag:
            new_source = create_new_Source(ra_target, dec_target, pos_err, flux, flux_err, band_fluxes, band_fluxerr,
                                           date)
            old_ms = list(load_specific_master_sources(session, tab_ms_id[index_impacted_master_sources[0]],obsid, ra_target, dec_target).values())[0]
            new_ms = MasterSource(session, old_ms.id,  list(old_ms.sources.values())+[new_source], old_ms.ra, old_ms.dec, old_ms.pos_err,[])
            new_ms.xmm_ul, new_ms.xmm_ul_dates = old_ms.xmm_ul,old_ms.xmm_ul_dates
            new_ms.slew_ul, new_ms.slew_ul_dates = old_ms.slew_ul,old_ms.slew_ul_dates
            new_ms.has_short_term_var= (old_ms.has_short_term_var or var_flag)
            # new_ms.simbad_type=old_ms.simbad_type[0]
            new_ms.simbad_type, new_ms.simbad_name = match_Simbad(ra_target, dec_target, pos_err)
            new_ms.var_ratio=var_ratio
            tab_alerts.append(new_ms)
            if var_ratio > tab_max[index_impacted_master_sources[0]]/mini_hist:
                if flux-flux_err > tab_max[index_impacted_master_sources[0]]:
                    flag_alerts.append('High Flux State')
                elif flux+flux_err < mini_hist:
                    flag_alerts.append('Low Flux State')
            else:
                flag_alerts.append('Past Variability')
            info_source.append(var_ratio)
            info_source.append(new_ms.simbad_type)
            info_source.append(new_ms.simbad_name)
    elif flag_known_ms!=-1: #We require the archival matching to not be ambiguous. If it was ambiguous, bad idea to compute the UpperLimits
        if ul:
            try:
                xmm_ul, xmm_ul_dates, slew_ul, slew_ul_dates = compute_upper_limit(ra_target, dec_target, flux, date)
                if len(xmm_ul+slew_ul)>0:
                    if (np.nanmin(xmm_ul+slew_ul) < (flux-flux_err)/5) or var_flag:
                        new_source = create_new_Source(ra_target, dec_target, pos_err, flux, flux_err, band_fluxes, band_fluxerr,
                                                       date)
                        new_ms = MasterSource(session, - 1, [new_source], new_source.ra, new_source.dec, new_source.poserr, [])
                        new_ms.has_short_term_var = var_flag
                        new_ms.xmm_ul = xmm_ul
                        new_ms.xmm_ul_dates = xmm_ul_dates
                        new_ms.slew_ul = slew_ul
                        new_ms.slew_ul_dates = slew_ul_dates
                        new_ms.var_ratio = (flux-flux_err)/np.nanmin(xmm_ul+slew_ul)
                        new_ms.simbad_type, new_ms.simbad_name = match_Simbad(ra_target,dec_target,pos_err)
                        tab_alerts.append(new_ms)
                        flag_alerts.append('First Detection')
                        info_source.append((flux-flux_err)/np.nanmin(xmm_ul+slew_ul))
                        info_source.append(new_ms.simbad_type)
                        info_source.append(new_ms.simbad_name)
            except:
                pass
        elif var_flag:
            new_source = create_new_Source(ra_target, dec_target, pos_err, flux, flux_err, band_fluxes, band_fluxerr,
                                           date)
            new_ms = MasterSource(session, - 1, [new_source], new_source.ra, new_source.dec, new_source.poserr, [])
            new_ms.has_short_term_var = var_flag
            new_ms.simbad_type = "Not Checked"
            new_ms.simbad_name = "Not Checked"
            tab_alerts.append(new_ms)
            flag_alerts.append('Short-term Variable Detection')
            info_source.append(0)
            info_source.append(new_ms.simbad_type)
            info_source.append(new_ms.simbad_name)


    return tab_alerts, flag_alerts, info_source


def testing_functions_NGC7793():
    """
    Testing function, meant to work on a sub-sample of the catalog corresponding to an observation of NGC 7793 containing
    several variables objects
    Might not work anymore because of the dictionary metadata update.
    """
    #Initial input to get the Ra, Dec, Flux, FluxErr
    pos_err = 1#input("1sigma PosErr of the source (in arcsec)?")
    flux_level = 1e-14#input("Flux Level of the source ?")
    flux_err_level = 0.5e-14#input("Flux Error Level of the source ?")
    band_fluxes = [[1e-12]*5]
    band_fluxerr = [[[3e-13]*5],[[3e-13]*5]]
    date = Time(2022.0, format="decimalyear").mjd


    raw_data = fits.open(f"{PATHTO.master_sources}Master_source_TestPipeline_NGC7793.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)

    tab_times = []
    tab_band_fluxes = [[list(line)] for line in sources_raw["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"]]
    tab_band_fluxerr = []
    for line in sources_raw["EP_1_FLUX_ERR","EP_2_FLUX_ERR","EP_3_FLUX_ERR","EP_4_FLUX_ERR","EP_5_FLUX_ERR"]:
        tab_band_fluxerr.append([[list(line)], [list(line)]])

    tab_times=[]
    tab_alerts=[]
    start = time.time()
    tab_alerts += transient_alert(1, 359.38, -32.584, 1, 2e-12, 1e-13, band_fluxes,
                                  band_fluxerr, date, var_flag=True)
    end = time.time()

    tab_times.append(end - start)
    pbar=tqdm(total=len(sources_raw))
    for ra, dec, flux, flux_err, band_flux, band_fluxerr, date in \
            zip(sources_raw["MS_RA"],sources_raw["MS_DEC"], sources_raw["EP_8_FLUX"],\
                    sources_raw["EP_8_FLUX_ERR"],tab_band_fluxes, tab_band_fluxerr,[56061.15969907407]*len(sources_raw)):
        start = time.time()
        tab_alerts += transient_alert(1, ra, dec, 5, flux, flux_err, band_flux,
                                     band_fluxerr, date, var_flag=False)[0]
        end = time.time()
        tab_times.append(end - start)
        pbar.update(1)
    pbar.close()

    plt.hist(tab_times, bins=np.geomspace(1e-3,1e2,20))
    plt.xscale("log")
    for ms in tab_alerts:
        ms.save_lightcurve({"ObsID":obsid})
#testing_functions_NGC7793()


ra=0 #in Degrees
dec=0 #in Degrees
pos_err = 1  # input("1sigma PosErr of the source (in arcsec)?")
flux_level = 1e-14  # input("Flux Level of the source ?")
flux_err_level = 0.5e-14  # input("Flux Error Level of the source ?")
band_fluxes = [[1e-12] * 5]
band_fluxerr = [[[3e-13] * 5], [[3e-13] * 5]]
date = Time(2022.0, format="decimalyear").mjd
var_flag=True
obsid=1
#result = transient_alert(obsid, ra, dec, pos_err, flux_level, flux_err_level, band_fluxes,band_fluxerr, date,var_flag=False, ul=True)
#for ms in result:
#    ms.save_lightcurve(obsid)
