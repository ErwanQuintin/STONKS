"""
The Search for Transient Objects in New detections using Known Sources (STONKS) software is meant to allow automatic
comparison between the PPS file of a new XMM-Newton observation, and all the available archival X-ray data. It will send
out an alert and plot the lightcurve and available data for all sources with a long-term variability over a factor of 5.

This requires the script to be run with the same file structure as on the GitHub, to ensure access to hand-made catalog
data. You also need to state the path to the PPS folder.

The LoadMasterSources.py module loads the archival X-ray MasterSource objects, which correspond to associations between
various catalog sources. It will then be compared to new detection of a given PPS file in the STONKS_pipeline_alert.py
module.

It is also possible to use the archival catalog strictly for data mining past data, but this is done in the other script
(StudyMasterSources.py).

This software was developed as part of the XMM2ATHENA project. This project has received funding from the European
Union's Horizon 2020 research and innovation programme under grant agreement n°101004168, the XMM2ATHENA project.

Author: Erwan Quintin, erwan.quintin@irap.omp.eu
"""


import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from tqdm import tqdm
from astropy.constants import c
from scipy.stats import pearsonr
import webbrowser

style="bmh"
cmap_to_use="turbo"

def click_action(ra, dec, xmm_name, swift_name):
    """
    Function that triggers when clicking on the lightcurve plot. Opens several browser pages at the position (ra,dec).
    The catalogs name are used, if != [], to open the corresponding catalog pages
    :param ra: RA of the source
    :param dec: DEC of the source
    :param xmm_name: IAUNAME of the XMM-Newton source, [] if no XMM-Newton counterpart
    :param swift_name: IAUNAME of the Swift source, [] if no Swift counterpart
    :return: Nothing is returned
    """
    url_esaskyDSS = "http://sky.esa.int/?target="+str(np.round(ra,4))+" "+str(np.round(dec,4))+"&hips=DSS2+color&fov=0.1&cooframe=J2000&sci=true&lang=en"
    url_esaskyXMM = "http://sky.esa.int/?target=" + str(np.round(ra, 4)) + " " + str(
        np.round(dec, 4)) + "&hips=XMM-Newton+EPIC+color&fov=0.1&cooframe=J2000&sci=true&lang=en"
    url_esaskyChandra = "http://sky.esa.int/?target=" + str(np.round(ra, 4)) + " " + str(
        np.round(dec, 4)) + "&hips=Chandra+RGB&fov=0.1&cooframe=J2000&sci=true&lang=en"
    webbrowser.open_new(url_esaskyDSS)
    webbrowser.open(url_esaskyXMM, new=0)
    webbrowser.open(url_esaskyChandra, new=0)
    if xmm_name != []:
        xmm_name = xmm_name[5:].replace(' ', '%2B')
        url_xmmssc = f"http://xmm-catalog.irap.omp.eu/sources?f={xmm_name}"
        webbrowser.open(url_xmmssc, new=0)
    if swift_name != []:
        url_swift = "https://www.swift.ac.uk/2SXPS/" + swift_name
        webbrowser.open(url_swift, new=0)
    url_simbad = "http://simbad.u-strasbg.fr/simbad/sim-coo?Coord="+str(ra)+"+"+str(dec)+"&Radius=1&Radius.unit=arcmin&submit=submit+query"
    webbrowser.open(url_simbad)
    url_rosat= f"http://xmm-ssc.irap.omp.eu/claxson/xray_analyzer2.php?srcquery={ra}%20{dec}"
    webbrowser.open(url_rosat)


#Starting from here, we define all the relevant column names or catalog properties that will be used later on.
#If you're not interested in a catalog for whatever reason, you can remove it. They need to stay in the overall same order.
#The NewXMM catalog is a placeholder, used for new detections for the pipeline alert system. It will have the same
#properties as XMM where needed.

catalogs = ["XMM","Chandra","Swift","eRosita","Slew","RASS","WGACAT","Stacked","NewXMM"]

posErr_Names = {}
src_names={}
colors = {}
for ind,cat in enumerate(catalogs):
    posErr_Names[cat]=f"{cat}_PosErr"
    src_names[cat] = f"{cat}_IAUNAME"
    colors[cat] = matplotlib.cm.get_cmap(cmap_to_use)(ind / len(catalogs)) #Used for lightcurve plotting


optical_catalogs = ["OM","UVOT"]
for opt_cat in optical_catalogs:
    src_names[opt_cat]=f"{opt_cat}_IAUNAME"
optical_formats={'OM':'o','UVOT':'s'}
optical_colors={}
for ind, band in enumerate(["UVW2","UVM2","UVW1","U","B","V"]):
    optical_colors[band] = matplotlib.cm.get_cmap(cmap_to_use)(ind / 6)
obsid_names={"XMM":"OBS_ID","Swift":"ObsID", "Stacked":"OBS_ID", "OM":"OBSID","UVOT":"OBSID"}
short_term_var_name={"XMM":"VAR_FLAG","Chandra":"var_code"}


flux_names={"XMM": "EP_8_FLUX",
            "Chandra":"flux_aper_b",
            "Swift":"Flux",
            "eRosita":"ML_FLUX",
            "Slew":"Flux",
            "Stacked":"EP_FLUX",
            "RASS":"Flux",
            "WGACAT":"Flux"}
flux_error_names={"XMM": ["EP_8_FLUX_ERR","EP_8_FLUX_ERR"],
                  "Chandra":["flux_aper_b_negerr","flux_aper_b_poserr"],
                  "Swift":["FluxErr_neg","FluxErr_pos"],
                  "eRosita":["ML_FLUX_ERR","ML_FLUX_ERR"],
                  "Slew":["FluxErr","FluxErr"],
                  "Stacked":["EP_FLUX_ERR","EP_FLUX_ERR"],
                  "RASS":["FluxErr","FluxErr"],
                  "WGACAT":["FluxErr","FluxErr"]}
conv_factors = {"XMM": 1/0.999,
                "NewXMM": 1/0.999,
                "Chandra":1/0.69,
                "Swift":1/0.9,
                "eRosita":1/0.39,
                "Slew":1/0.999,
                "Stacked":1/0.999,
                "RASS":1/0.35,
                "WGACAT":1/0.35}
time_names={"XMM": "MJD_START",
            "Chandra":"gti_mjd_obs",
            "Swift":"MidTime_MJD",
            "eRosita":"MJD_OBS",
            "Slew":"DATE_OBS",
            "Stacked":"MJD_FIRST",
            "RASS":"OBS_DATE_1",
            "WGACAT":"StartDate",
            "OM":"MJD_START",
            "UVOT":"DATE_MIN"}

band_flux_names = {"XMM":["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"],
                  "Chandra":["flux_aper_s","flux_aper_m","flux_aper_h"],
                  "Swift":["Flux1","Flux2","Flux3"],
                  "Slew":["Flux6","Flux7"],
                  "eRosita":['ML_FLUX_b1','ML_FLUX_b2','ML_FLUX_b3','ML_FLUX_b4'],
                  "RASS":["Flux1","Flux3","Flux4"],
                  "WGACAT":["Flux1","Flux2","Flux3"],
                  "Stacked":["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"],
                  "OM":["UVW2_AB_FLUX","UVM2_AB_FLUX","UVW1_AB_FLUX","U_AB_FLUX","B_AB_FLUX","V_AB_FLUX"],
                  "UVOT":["UVW2_FLUX","UVM2_FLUX","UVW1_FLUX","U_FLUX","B_FLUX","V_FLUX"]}
band_fluxerr_names = {"XMM":[["EP_1_FLUX_ERR","EP_2_FLUX_ERR","EP_3_FLUX_ERR","EP_4_FLUX_ERR","EP_5_FLUX_ERR"],
                            ["EP_1_FLUX_ERR","EP_2_FLUX_ERR","EP_3_FLUX_ERR","EP_4_FLUX_ERR","EP_5_FLUX_ERR"]],
                     "Chandra":[["flux_aper_s_negerr","flux_aper_m_negerr","flux_aper_h_negerr"],
                                ["flux_aper_s_poserr","flux_aper_m_poserr","flux_aper_h_poserr"]],
                     "Swift":[["FluxErr1_neg","FluxErr2_neg","FluxErr3_neg"],["FluxErr1_pos","FluxErr2_pos","FluxErr3_pos"]],
                     "Slew":[["Flux6Err","Flux7Err"],["Flux6Err","Flux7Err"]],
                     "eRosita":[['ML_FLUX_ERR_b1','ML_FLUX_ERR_b2','ML_FLUX_ERR_b3','ML_FLUX_ERR_b4'],
                                      ['ML_FLUX_ERR_b1','ML_FLUX_ERR_b2','ML_FLUX_ERR_b3','ML_FLUX_ERR_b4']],
                     "RASS":[["FluxErr1","FluxErr3","FluxErr4"],["FluxErr1","FluxErr3","FluxErr4"]],
                     "WGACAT":[["FluxErr1","FluxErr2","FluxErr3"],["FluxErr1","FluxErr2","FluxErr3"]],
                     "Stacked": [["EP_1_FLUX_ERR", "EP_2_FLUX_ERR", "EP_3_FLUX_ERR", "EP_4_FLUX_ERR", "EP_5_FLUX_ERR"],
                              ["EP_1_FLUX_ERR", "EP_2_FLUX_ERR", "EP_3_FLUX_ERR", "EP_4_FLUX_ERR", "EP_5_FLUX_ERR"]],
                     "OM":["UVW2_AB_FLUX_ERR","UVM2_AB_FLUX_ERR","UVW1_AB_FLUX_ERR","U_AB_FLUX_ERR","B_AB_FLUX_ERR","V_AB_FLUX_ERR"],
                     "UVOT":["UVW2_FLUX_ERR","UVM2_FLUX_ERR","UVW1_FLUX_ERR","U_FLUX_ERR","B_FLUX_ERR","V_FLUX_ERR"]
                      }
band_center = {"XMM":[0.35,0.75,1.5,3.25,8.25],
               "NewXMM":[0.35,0.75,1.5,3.25,8.25],
               "Chandra":[0.85,1.6,4.5],
               "Swift":[0.65,1.5,6],
               "Slew":[1.1,7],
               "eRosita":[0.35,0.75,1.5,3.25],
               "RASS":[0.25,0.7,1.45],
               "WGACAT":[0.25,0.7,1.45],
               "Stacked":[0.35,0.75,1.5,3.25,8.25],
               "OM":[2000,2300,2750,3500,4500,5500],
               "UVOT":[2000,2300,2750,3500,4500,5500]}
band_half_width = {"XMM":[0.15,0.25,0.5,1.25,3.75],
              "NewXMM":[0.15,0.25,0.5,1.25,3.75],
              "Chandra":[0.35,0.4,2.5],
              "Swift":[0.35,0.5,4],
              "Slew":[0.9,5],
              "eRosita":[0.15,0.25,0.5,1.25],
              "RASS":[0.15,0.2,0.55],
              "WGACAT":[0.15,0.2,0.55],
              "Stacked":[0.15,0.25,0.5,1.25,3.75],
              "OM":[100,200,250,500,500,500],
              "UVOT":[100,200,250,500,500,500]}

optical_effective_wavelengths=[2120,2310,2910,3440,4500,5430] #In Angstroms, used to convert from erg/s/cm2/angstroms to erg/s/cm2



band_edges = {}
frequencies = {}
frequencies_half_width = {}
for cat in band_center.keys():
    band_edges[cat] = [center-width for (center, width) in zip(band_center[cat],band_half_width[cat])]
    band_edges[cat].append(band_center[cat][-1]+band_half_width[cat][-1])
for cat in catalogs:
    frequencies[cat]=[2.41e17*center for center in band_center[cat]]
    frequencies_half_width[cat] = [2.41e17*width for width in band_half_width[cat]]
xband_average_frequency = 2 * 2.41e17 #2keV in Hz, to compute alpha_OX
xband_width = 11.9 * 2.41e17 #11.9 keV in Hz, to compute alpha_OX
for cat in optical_catalogs:
    frequencies[cat] = [c.to('angstrom/s').value/wavelength for wavelength in optical_effective_wavelengths]
    frequences_max = [c.to('angstrom/s').value/wavelength for wavelength in [1800,1990,2430,3030,3810,5020]]
    frequences_min = [c.to('angstrom/s').value/wavelength for wavelength in [2550,2700,3610,3890,4900,5870]]
    frequencies_half_width[cat]  = [(freq_max-freq_min)/2 for (freq_max, freq_min) in zip(frequences_max,frequences_min)]

#Allows to define the limiting index between what is considered soft band (<2keV) and hard band, for each instrument
hr_bandlimit_index = {"XMM":3,"NewXMM":3,"Chandra":2,"Swift":2,"Slew":1,"eRosita":3,"RASS":3,"WGACAT":3,"Stacked":3}

#Conversion factors used to convert the flux in the instrument full band to the 0.1-12 keV band, assuming a spectral
#shape given by a standard powerlaw of photon index Gamma=1.7 and absorption nH=3e20 cm-2.
band_conv_factors_soft = {"XMM":0.35/0.35,
                          "NewXMM":0.35/0.35,
                          "Chandra":0.35/0.28,
                          "Swift":0.35/0.34,
                          "Slew":0.35/0.35,
                          "eRosita":0.35/0.35,
                          "RASS":0.35/0.35,
                          "WGACAT":0.35/0.35,
                          "Stacked":0.35/0.35}
band_conv_factors_hard = {"XMM":0.65/0.65,
                          "NewXMM":0.65/0.65,
                          "Chandra":0.65/0.41,
                          "Swift":0.65/0.56,
                          "Slew":0.65/0.65,
                          "eRosita":0.65/0.25,
                          "RASS":np.nan,
                          "WGACAT":np.nan,
                          "Stacked":0.65/0.65}

hr_track_markers = {"XMM":"o","NewXMM":'o',"Chandra":"v","Swift":"s","Slew":"P","eRosita":"^","RASS":"d","WGACAT":"d","Stacked":"*"}

path_to_catalogs = "Data/Catalogs/FullData/"
path_to_master_sources = "Data/MasterSource/"

class Source:
    """
    A Source object corresponds to a source from one of the X-ray catalogs. It has several attributes:
    - catalog: the corresponding catalog name, in the same naming convention as the catalog Table defined at the top
    - iau_name: the name of the source, considered as a unique identifier
    - fluxes: a table containing the fluxes extrapolated in the 0.1-12keV band
    - flux_errors: a table containing the 1 sigma flux errors extrapolated in the 0.1-12keV band. It consists of 2 tables,
    one for negative errors and one for positive errors.
    - timesteps: a table containing the MJD dates of detections
    - obsids: a table containing the ObsID for each detection, in the case of XMM and Swift (used for matching with OM & UVOT)
    - band_flux and band_fluxerr: same as fluxes and flux_errors, but is divided in the various detection bands of each instrument.
    Used for spectrum & SED plotting

    In the end, each Source will be associated to a unique MasterSource, each MasterSource having Source objects from several distinct catalogs
    """

    def __init__(self, catalog, iau_name, flux, fluxerr, timesteps, band_flux, band_fluxerr, obsids=[],swift_stacked_flux=[],
                 swift_stacked_flux_err=[[],[]], swift_stacked_times=[[],[]], xmm_offaxis=[], short_term_var=[]):
        """
        Initialisation function, used to build the Source object. It will also compute Hardness Ratios using the fluxes
        in soft and hard bands, if available, as well as the error on this HR.
        For Swift, we have the stacked detection, which is treated separately, and corresponds to the average flux over
        all observations, even in non detections; it is a proxy to upper limits if it is constraining.
        For XMM-Newton, we add information about the off-axis angle and the short term variability.
        A variability amplitude is also computed for the Source.
        """
        self.catalog = catalog
        self.name = iau_name
        self.master_source = []
        self.fluxes = flux
        self.flux_errors=fluxerr
        self.timesteps=[float(elt) for elt in timesteps]
        self.obsids=[int(obsid) for obsid in obsids]

        self.band_flux = band_flux
        self.band_fluxerr = band_fluxerr

        self.soft_dets = [np.sum(det[:hr_bandlimit_index[catalog]])*band_conv_factors_soft[catalog] for det in self.band_flux]
        self.soft_errors = [[np.sqrt(np.sum(np.array(err_neg[:hr_bandlimit_index[catalog]])**2))*band_conv_factors_soft[catalog] for err_neg in self.band_fluxerr[0]],
                            [np.sqrt(np.sum(np.array(err_pos[:hr_bandlimit_index[catalog]])**2))*band_conv_factors_soft[catalog] for err_pos in self.band_fluxerr[1]]]
        if catalog!= "RASS" and catalog!="WGACAT":
            self.hard_dets = [np.sum(det[hr_bandlimit_index[catalog]:])*band_conv_factors_hard[catalog] for det in self.band_flux]
            self.hard_errors = [
                [np.sqrt(np.sum(np.array(err_neg[hr_bandlimit_index[catalog]:])**2)) * band_conv_factors_hard[catalog] for err_neg in
                 self.band_fluxerr[0]],
                [np.sqrt(np.sum(np.array(err_pos[hr_bandlimit_index[catalog]:])**2)) * band_conv_factors_hard[catalog] for err_pos in
                 self.band_fluxerr[1]]]
        else:
            self.hard_dets = [np.nan for det in self.fluxes]
            self.hard_errors = [[np.nan for det in self.fluxes],[np.nan for det in self.fluxes]]


        self.hardness = [(hard-soft)/(hard+soft) for (soft,hard) in zip(self.soft_dets, self.hard_dets)]
        low_soft = np.where(np.array(self.soft_dets) - np.array(self.soft_errors[0]) < 0, 0,
                            np.array(self.soft_dets) - np.array(self.soft_errors[0]))
        low_hard = np.where(np.array(self.hard_dets) - np.array(self.hard_errors[0]) < 0, 0,
                            np.array(self.hard_dets) - np.array(self.hard_errors[0]))
        up_soft = np.where(np.array(self.soft_dets) + np.array(self.soft_errors[1]) < 0, 0,
                           np.array(self.soft_dets) + np.array(self.soft_errors[1]))
        up_hard = np.where(np.array(self.hard_dets) + np.array(self.hard_errors[1]) < 0, 0,
                           np.array(self.hard_dets) + np.array(self.hard_errors[1]))
        self.hardness_err = [[hr - (hard-soft)/(hard+soft) for (soft,hard,hr) in zip(up_soft, low_hard, self.hardness)],
                             [(hard-soft)/(hard+soft) - hr for (soft,hard,hr) in zip(low_soft, up_hard, self.hardness)]]
        self.swift_stacked_flux = swift_stacked_flux
        self.swift_stacked_flux_err = swift_stacked_flux_err
        self.swift_stacked_times = swift_stacked_times
        self.swift_stacked_variable = False

        self.min_upper = 1
        self.max_lower = 0
        self.var = 1
        if len(flux)>0:
            self.min_upper = min(np.array(flux) + np.array(fluxerr[1]))
            self.max_lower = max(np.array(flux) - np.array(fluxerr[0]))
        if swift_stacked_flux!=[]:
            stacked_min = min(np.array(swift_stacked_flux)+np.array(swift_stacked_flux_err[1]))
            if stacked_min<0.5*self.min_upper:
                self.swift_stacked_variable = True
            self.min_upper = min(self.min_upper, stacked_min)
        if len(flux)+len(swift_stacked_flux) > 1:
            self.var = self.max_lower/self.min_upper

        self.xmm_offaxis = xmm_offaxis
        self.short_term_var = short_term_var

class OpticalSource:
    """
    An OpticalSource is the equivalent of a Source, but for the OM & UVOT catalogs. It has less properties. Here, the
    bands are the 6 optical filters. The ObsIDs are used for matching with X-rays. In the end, an OpticalSource is
    associated to a MasterSource, which contains also the X-ray sources. It is possible with this framework to work
    solely on the optical catalogs, but this work will be done later.
    """

    def __init__(self, catalog, iau_name, timesteps, band_flux, band_fluxerr, obsids):
        """
        Initialisation function, used to create an OpticalSource object.
        """
        self.catalog = catalog
        self.name = iau_name
        self.master_source = []
        self.band_flux = band_flux
        self.band_fluxerr=band_fluxerr
        self.timesteps=[float(elt) for elt in timesteps]
        self.obsids = [int(obsid) for obsid in obsids]

class MasterSource:
    """
    A MasterSource corresponds to a single physical source, built on the association of multiple archival catalogs.
    A MasterSource has several attributes:
    - source: A dictionary which gives access to the underlying catalogs sources, which are Source objects in our framework.
    The keys of this dictionary are the names of the corresponding catalogs.
    - source_fluxes and source_error_bars: compiled flux and flux_errors from all its constituting Source objects, it is
    useful in order to access overall flux properties of the MasterSource (maximum flux, variability,...).
    - tab_hr and tab_hr_err: same thing but for hardness properties instead of flux properties.
    - var_ratio, var_amplitude, var_significance: correspond to different ways of quantifying flux variability of the source
    - hr_var, hr_var_signif: same thing but for hardness properties instead of flux properties

    A MasterSource only has one method, plot_lightcurve(), which produces a multi-panel plot of all relevant information
    """
    def __init__(self, id, tab_sources, ra, dec, poserr, tab_optical_sources):
        """
        Initialisation function, used to build a MasterSource object. We also compile the multi-instrument properties at
        this stage (variabilities, HR variabilities,...)
        :param id: Identifier of the MasterSource, used to access it in a dictionary with ms.id as a key, and ms as value
        :param tab_sources: Table containing all the catalog Source objects
        :param ra: RA of the MasterSource computed as weighted average of the constituting Source objects
        :param dec: Dec of the MasterSource computed as weighted average of the constituting Source objects
        :param poserr: 1 sigma Position Error of the MasterSource computed as weighted average of the constituting Source objects
        :param tab_optical_sources: Table containing the OpticalSource objects, if any
        """
        self.id = id
        self.sources = {}
        self.sources_fluxes = []
        self.sources_error_bars = []
        self.sources_var = []
        self.tab_hr = []
        self.tab_hr_err = [[],[]]
        self.never_on_axis_xmm = False
        self.has_short_term_var = False
        self.min_time=60000
        self.max_time=0
        for source in tab_sources:
            if ("XMM" in self.sources.keys()) and (source.catalog == "Stacked"):
                #We remove the Stacked detection that correspond to a clean XMM detection
                xmm_obsid = self.sources["XMM"].obsids
                stacked_obsid = source.obsids
                new_det_ind = [i for i in range(len(stacked_obsid)) if stacked_obsid[i] not in xmm_obsid]
                source.fluxes = source.fluxes[new_det_ind]
                source.flux_errors[0] = source.flux_errors[0][new_det_ind]
                source.flux_errors[1] = source.flux_errors[1][new_det_ind]
                source.timesteps = np.array(source.timesteps)[new_det_ind]
                source.obsids = np.array(source.obsids)[new_det_ind]
                source.hardness = np.array(source.hardness)[new_det_ind]
                source.hardness_err[0] = np.array(source.hardness_err[0])[new_det_ind]
                source.hardness_err[1] = np.array(source.hardness_err[1])[new_det_ind]

                source.band_flux = source.band_flux[new_det_ind]
                source.band_fluxerr[0] = source.band_fluxerr[0][new_det_ind]
                source.band_fluxerr[1] = source.band_fluxerr[1][new_det_ind]
            source.master_source = self
            self.sources[source.catalog]=source
            for (flux, fluxerr_neg, fluxerr_pos) in zip(source.fluxes, source.flux_errors[0], source.flux_errors[1]):
                self.sources_fluxes.append(flux)
                self.sources_error_bars.append(max(fluxerr_neg, fluxerr_pos))
                self.sources_var.append(source.var)
            self.tab_hr += list(source.hardness)
            self.tab_hr_err[0] += list(source.hardness_err[0])
            self.tab_hr_err[1] += list(source.hardness_err[1])
            for (flux, fluxerr_neg, fluxerr_pos, start, stop) in zip(source.swift_stacked_flux, source.swift_stacked_flux_err[0], source.swift_stacked_flux_err[1], source.swift_stacked_times[0], source.swift_stacked_times[1]):
                self.sources_fluxes.append(flux)
                self.sources_error_bars.append(max(fluxerr_neg, fluxerr_pos))
                self.min_time = min(start, self.min_time)
                self.max_time = max(stop, self.max_time)
            if source.xmm_offaxis!=[]:
                if np.nanmin(source.xmm_offaxis)>1:
                    self.never_on_axis_xmm = True
            if source.timesteps!=[]:
                self.min_time = min(min(source.timesteps), self.min_time)
                self.max_time = max(max(source.timesteps), self.max_time)
            for var_flag in source.short_term_var:
                if var_flag>0:
                    self.has_short_term_var=True
        self.sources_fluxes = np.array(self.sources_fluxes)
        self.sources_error_bars = np.array(self.sources_error_bars)


        #Min upper and Max lower correspond respectively to the minimum of the fluxes with added error bars, and to the
        #maximum of the fluxes with subtracted error bars. It allows for a conservative estimate of variability,
        #that can thus be below one in case it is consistent with being constant within error bars.
        self.min_upper = 1
        self.max_lower = 0
        self.var_ratio = 1
        self.var_amplitude = 0
        self.var_significance = 0
        if len(self.sources_fluxes)>0 and (not np.isnan(self.sources_fluxes).all()):
            min_upper_ind = np.argmin(self.sources_fluxes + self.sources_error_bars)
            self.min_upper = (self.sources_fluxes + self.sources_error_bars)[min_upper_ind]
            max_lower_ind = np.argmax(self.sources_fluxes - self.sources_error_bars)
            self.max_lower = (self.sources_fluxes - self.sources_error_bars)[max_lower_ind]
            self.var_ratio = self.max_lower/self.min_upper
            self.var_amplitude = self.max_lower - self.min_upper
            self.var_optimistic = self.sources_fluxes[max_lower_ind]/self.sources_fluxes[min_upper_ind]
            self.var_significance = self.var_amplitude/np.sqrt(self.sources_error_bars[max_lower_ind]**2 + self.sources_error_bars[min_upper_ind]**2)

        self.hr_min = np.nan
        self.hr_max = np.nan
        self.hr_var = np.nan
        self.hr_var_signif = np.nan
        if len(self.tab_hr)>1 and (not np.isnan(self.tab_hr).all()) and (not np.isnan(self.tab_hr_err).all()):
            index_hr_min = np.nanargmin(np.array(self.tab_hr)+np.array(self.tab_hr_err[1]))
            index_hr_max = np.nanargmax(np.array(self.tab_hr)-np.array(self.tab_hr_err[0]))
            self.hr_min = (np.array(self.tab_hr)+np.array(self.tab_hr_err[1]))[index_hr_min]
            self.hr_max = (np.array(self.tab_hr)-np.array(self.tab_hr_err[0]))[index_hr_max]
            self.hr_var = self.hr_max - self.hr_min
            if self.tab_hr_err[1][index_hr_min]**2 + self.tab_hr_err[0][index_hr_max]**2 >0:
                self.hr_var_signif = self.hr_var/np.sqrt(self.tab_hr_err[1][index_hr_min]**2 + self.tab_hr_err[0][index_hr_max]**2)
            else:
                self.hr_var_signif = np.nan

        #Upper limits from XMM, Slew, and Chandra, will be loaded later on
        self.xmm_ul = []
        self.xmm_ul_dates = []
        self.xmm_ul_obsids = []

        self.slew_ul = []
        self.slew_ul_dates = []
        self.slew_ul_obsids = []

        self.chandra_ul = []
        self.chandra_ul_dates = []

        self.ra = float(ra)
        self.dec = float(dec)
        self.pos_err = float(poserr)

        #The variabilities of the OM & UVOT combined counterparts must be computed independently in each of the 6 filters
        self.optical_sources={}
        self.optical_min_upper = [1,1,1,1,1,1]
        self.optical_var = [1, 1, 1, 1, 1, 1]
        self.optical_max_lower = [0,0,0,0,0,0]
        for opt_source in tab_optical_sources:
            self.optical_sources[opt_source.catalog]=opt_source
            opt_source.master_source = self
            if len(opt_source.band_flux)>1:
                lower_fluxes = opt_source.band_flux-opt_source.band_fluxerr
                upper_fluxes = opt_source.band_flux+opt_source.band_fluxerr
                self.optical_min_upper = np.nanmin([self.optical_min_upper, np.nanmin(upper_fluxes, axis=0)], axis=0)
                self.optical_max_lower = np.nanmax([self.optical_max_lower, np.nanmax(lower_fluxes, axis=0)], axis=0)
                self.optical_var = self.optical_max_lower/self.optical_min_upper
            self.min_time = min(min(opt_source.timesteps), self.min_time)
            self.max_time = max(max(opt_source.timesteps), self.max_time)

        self.glade_distance=[]

        self.simbad_type=''

    def plot_lightcurve(self, limit_dates=[]):
        """
        Produces a multi-panel plot with most of the useful multi-instrument information about the source. From left to
        right and top to bottom:
        1. Long term multi-instrument X-ray lightcurves, extrapolated to the 0.1-12 keV
        2. Multi-instrument X-ray spectra, used to assess if the apparent flux change are due to a wrong extrapolation
        (if the spectra are comparable in common bands, but extrapolated fluxes are not)
        3. Hardness - Flux diagram
        4. Long term multi-instrument multi-band optical lightcurve from OM & UVOT data
        5. Multi-instrument SED, i.e. X-ray spectra and optical spectra associated
        6. Lightcurve of the slope of the SED, alpha_OX, in the different bands.
        :param limit_dates:
        :return:
        """
        plt.rcParams.update({'font.size': 8})
        if self.optical_sources!={}:
            fig, [[ax1, ax2, ax3], [ax4, ax5, ax6]] = plt.subplots(2, 3, figsize=(15,10))
        else:
            fig, [ax1, ax2, ax3] = plt.subplots(1,3, figsize=(15,5))

        plt.suptitle(f'Simbad type: {self.simbad_type}   -   More details', picker=True, bbox=dict(facecolor=(180 / 256., 204 / 256., 252 / 256.)))
        xmm_name=[]
        swift_name=[]
        if len(self.xmm_ul)!= 0:
            ax1.errorbar(self.xmm_ul_dates, self.xmm_ul, yerr=0.2 * np.array(self.xmm_ul), uplims=True, fmt='none', c=colors["XMM"], label="XMM non-det.")
        if len(self.slew_ul)!= 0:
            ax1.errorbar(self.slew_ul_dates, self.slew_ul, yerr=0.2 * np.array(self.slew_ul), uplims=True, fmt='none', c=colors["Slew"], label="Slew non-det.")
        if len(self.chandra_ul)!= 0:
            ax1.errorbar(self.chandra_ul_dates, self.chandra_ul, yerr=0.2 * np.array(self.chandra_ul), uplims=True, fmt='none', c=colors["Chandra"])


        hardness_track=[]
        hardness_err_track=[[],[]]
        luminosity_track=[]
        luminosity_err_track = [[], []]
        time_track=[]
        catalogs_track=[]
        for cat in catalogs:
            if cat in self.sources.keys():
                source = self.sources[cat]
                #Plot the X-ray full band lightcurves
                ax1.errorbar(np.array(source.timesteps), np.array(source.fluxes),
                            yerr=np.array(source.flux_errors), fmt="o",c=colors[cat],
                            label=source.name, markeredgecolor='gray')
                if cat == "Swift":
                    times=[(stop+start)/2 for (start,stop) in zip(source.swift_stacked_times[0],source.swift_stacked_times[1])]
                    timerange=[(stop-start)/2 for (start,stop) in zip(source.swift_stacked_times[0],source.swift_stacked_times[1])]
                    ax1.errorbar(times, source.swift_stacked_flux,
                                yerr=source.swift_stacked_flux_err, xerr=timerange,
                                 fmt="o", markeredgecolor='gray', c=colors[cat])
                tab_width = 2*np.array(band_half_width[cat])
                for det in range(len(source.band_flux)):
                    # Plot the X-ray spectra of each detection, using the fluxes in various bands
                    ax2.step(band_edges[cat], [source.band_flux[det][0]/tab_width[0]] + list(source.band_flux[det]/tab_width), c=colors[cat], where='pre')
                    ax2.errorbar(band_center[cat], source.band_flux[det]/tab_width,
                                 yerr=[source.band_fluxerr[0][det]/tab_width,source.band_fluxerr[1][det]/tab_width],
                                 fmt="o", markeredgecolor='gray', c=colors[cat], alpha=0.4)
                hardness_track+=list(source.hardness)
                hardness_err_track[0] += list(source.hardness_err[0])
                hardness_err_track[1] += list(source.hardness_err[1])
                luminosity_track+=list(source.fluxes)
                luminosity_err_track[0] += list(source.flux_errors[0])
                luminosity_err_track[1] += list(source.flux_errors[1])
                time_track+=list(source.timesteps)
                catalogs_track+=[cat for elt in source.timesteps]

                #Names are used for the web browser click button
                if cat=="XMM":
                    xmm_name = source.name
                elif cat=="Chandra":
                    chandra_name = source.name
                elif cat=="Swift":
                    swift_name = source.name

        #2nd line of the multi-panel plot is added only if there is an optical counterpart
        if self.optical_sources!={}:
            optical_band_observed={"UVW2":False,"UVM2":False,"UVW1":False,"U":False,"B":False,"V":False}
            for cat in optical_catalogs:
                if cat in self.optical_sources.keys():
                    opt_source = self.optical_sources[cat]
                    lightcurves = np.transpose(opt_source.band_flux)
                    ligthcurve_errors = np.transpose(opt_source.band_fluxerr)
                    for lightcurve, lightcurve_err, band in zip(lightcurves, ligthcurve_errors, ["UVW2","UVM2","UVW1","U","B","V"]):
                        if not np.isnan(lightcurve).all():
                            optical_band_observed[band]=True
                        #Plot the lightcurves in every filter
                        ax4.errorbar(opt_source.timesteps, lightcurve, yerr=lightcurve_err, fmt=optical_formats[cat], markeredgecolor='gray', c=optical_colors[band])

            for band in ["UVW2","UVM2","UVW1","U","B","V"]:
                if optical_band_observed[band]:
                    ax4.errorbar([], [], fmt="o", c=optical_colors[band], label=band, markeredgecolor='gray')

            #SED plot part
            all_freq = []
            all_nuFnu = []
            all_freq_width = []
            all_nuFnu_err = []
            all_formats = []
            all_names = []
            all_colors = []

            max_flux = np.nanmax(self.sources_fluxes)
            min_flux = np.nanmin(self.sources_fluxes)
            cmap = matplotlib.cm.get_cmap("RdYlBu_r")
            for (opt_cat, x_cat) in zip(["OM", "UVOT"], ['XMM', 'Swift']):
                if (opt_cat in self.optical_sources.keys()) and (x_cat in self.sources.keys()):
                    source = self.sources[x_cat]
                    x_obsid = source.obsids
                    opt_source = self.optical_sources[opt_cat]
                    opt_obsid = opt_source.obsids
                    common_obsids = set(x_obsid).intersection(opt_obsid)
                    indices_x = [x_obsid.index(obsid) for obsid in common_obsids]
                    indices_opt = [opt_obsid.index(obsid) for obsid in common_obsids]
                    indices_only_x = [x_obsid.index(obsid) for obsid in x_obsid if obsid not in common_obsids]
                    indices_only_opt = [opt_obsid.index(obsid) for obsid in opt_obsid if obsid not in common_obsids]
                    for ind_x, ind_opt in zip(indices_x, indices_opt):
                        #X-ray detections that have a matching optical detection
                        fluxes = list(opt_source.band_flux[ind_opt])[::-1] + list(elt for elt in source.band_flux[ind_x])
                        flux_errors = [list(opt_source.band_fluxerr[ind_opt])[::-1] + list(source.band_fluxerr[0][ind_x]),
                                       list(opt_source.band_fluxerr[ind_opt])[::-1] + list(source.band_fluxerr[1][ind_x])]

                        #To compute the SED we convert the energy into frequencies, and use these to have a nuFnu plot
                        tab_freq = frequencies[opt_cat][::-1] + frequencies[x_cat]
                        tab_freq_width = frequencies_half_width[opt_cat][::-1] + frequencies_half_width[x_cat]
                        nuFnu = [freq * flux / (2 * width) for (flux, freq, width) in zip(fluxes, tab_freq, tab_freq_width)]
                        nuFnu_err = [[freq * fluxerrneg / (2 * width) for (fluxerrneg, freq, width) in
                                      zip(flux_errors[0], tab_freq, tab_freq_width)],
                                     [freq * fluxerrpos / (2 * width) for (fluxerrpos, freq, width) in
                                      zip(flux_errors[1], tab_freq, tab_freq_width)]]
                        all_freq.append(tab_freq)
                        all_freq_width.append(tab_freq_width)
                        all_nuFnu.append(nuFnu)
                        all_nuFnu_err.append(nuFnu_err)
                        all_formats.append(hr_track_markers[x_cat])
                        all_names.append(source.name)
                        all_colors.append(
                            cmap(0.9 - np.log10(source.fluxes[ind_x] / min_flux) / np.log10(max_flux / min_flux)))
                    if len(indices_only_x) > 0:
                        #X-ray detections without optical counterpart
                        for ind_x in indices_only_x:
                            fluxes = list(elt for elt in source.band_flux[ind_x])
                            flux_errors = [list(source.band_fluxerr[0][ind_x]), list(source.band_fluxerr[1][ind_x])]
                            tab_freq = frequencies[x_cat]
                            tab_freq_width = frequencies_half_width[x_cat]
                            nuFnu = [freq * flux / (2 * width) for (flux, freq, width) in
                                     zip(fluxes, tab_freq, tab_freq_width)]
                            nuFnu_err = [[freq * fluxerrneg / (2 * width) for (fluxerrneg, freq, width) in
                                          zip(flux_errors[0], tab_freq, tab_freq_width)],
                                         [freq * fluxerrpos / (2 * width) for (fluxerrpos, freq, width) in
                                          zip(flux_errors[1], tab_freq, tab_freq_width)]]
                            all_freq.append(tab_freq)
                            all_freq_width.append(tab_freq_width)
                            all_nuFnu.append(nuFnu)
                            all_nuFnu_err.append(nuFnu_err)
                            all_formats.append(hr_track_markers[x_cat])
                            all_names.append(source.name)
                            all_colors.append(
                                cmap(0.9 - np.log10(source.fluxes[ind_x] / min_flux) / np.log10(max_flux / min_flux)))
                    if len(indices_only_opt) > 0:
                        #Optical detections without X-ray counterparts
                        for ind_opt in indices_only_opt:
                            fluxes = list(opt_source.band_flux[ind_opt])[::-1]
                            flux_errors = [list(opt_source.band_fluxerr[ind_opt])[::-1],
                                           list(opt_source.band_fluxerr[ind_opt])[::-1]]
                            tab_freq = frequencies[opt_cat][::-1]
                            tab_freq_width = frequencies_half_width[opt_cat][::-1]
                            nuFnu = [freq * flux / (2 * width) for (flux, freq, width) in
                                     zip(fluxes, tab_freq, tab_freq_width)]
                            nuFnu_err = [[freq * fluxerrneg / (2 * width) for (fluxerrneg, freq, width) in
                                          zip(flux_errors[0], tab_freq, tab_freq_width)],
                                         [freq * fluxerrpos / (2 * width) for (fluxerrpos, freq, width) in
                                          zip(flux_errors[1], tab_freq, tab_freq_width)]]
                            all_freq.append(tab_freq)
                            all_freq_width.append(tab_freq_width)
                            all_nuFnu.append(nuFnu)
                            all_nuFnu_err.append(nuFnu_err)
                            all_formats.append(hr_track_markers[x_cat])
                            all_names.append(source.name)
                            all_colors.append(cmap(0.5))

            for x_cat in catalogs:
                #We add the X-ray detections of catalogs other than XMM and Swift
                if (x_cat in self.sources.keys()) and (x_cat not in ['XMM', 'Swift']):
                    source = self.sources[x_cat]
                    for ind_x in range(len(source.band_flux)):
                        fluxes = source.band_flux[ind_x]
                        flux_errors = [source.band_fluxerr[0][ind_x], source.band_fluxerr[1][ind_x]]
                        tab_freq = frequencies[x_cat]
                        tab_freq_width = frequencies_half_width[x_cat]
                        nuFnu = [freq * flux / (2 * width) for (flux, freq, width) in zip(fluxes, tab_freq, tab_freq_width)]
                        nuFnu_err = [[freq * fluxerrneg / (2 * width) for (fluxerrneg, freq, width) in
                                      zip(flux_errors[0], tab_freq, tab_freq_width)],
                                     [freq * fluxerrpos / (2 * width) for (fluxerrpos, freq, width) in
                                      zip(flux_errors[1], tab_freq, tab_freq_width)]]
                        all_freq.append(tab_freq)
                        all_freq_width.append(tab_freq_width)
                        all_nuFnu.append(nuFnu)
                        all_nuFnu_err.append(nuFnu_err)
                        all_formats.append(hr_track_markers[x_cat])
                        all_names.append(source.name)
                        all_colors.append(
                            cmap(0.9 - np.log10(source.fluxes[ind_x] / min_flux) / np.log10(max_flux / min_flux)))

            #We plot the final obtained SEDs, adding a legend entry only for the first detection of each source
            seen_names = []
            for (ind, tab_freq, tab_freq_width, nuFnu, nuFnu_err, format, name, color) in zip(range(len(all_freq)),
                                                                                              all_freq, all_freq_width,
                                                                                              all_nuFnu, all_nuFnu_err,
                                                                                              all_formats, all_names,
                                                                                              all_colors):
                ax5.plot(np.array(tab_freq)[~np.isnan(nuFnu)], np.array(nuFnu)[~np.isnan(nuFnu)], c=color, ls="--",
                         alpha=0.4)
                if name not in seen_names:
                    ax5.errorbar(tab_freq, nuFnu, xerr=tab_freq_width, yerr=nuFnu_err, fmt=format, label=name, c=color)
                    seen_names.append(name)
                else:
                    ax5.errorbar(tab_freq, nuFnu, xerr=tab_freq_width, yerr=nuFnu_err, fmt=format, c=color)

            #We plot the Alpha_OX lightcurves
            for band in ["UVW2","UVM2","UVW1","U","B","V"]:
                if self.alpha_band_x_Fnu[band] != []:
                    ax6.errorbar(self.alpha_band_x_timesteps[band], self.alpha_band_x_Fnu[band], yerr=np.transpose(self.alpha_band_x_Fnu_error[band]), fmt="o", markeredgecolor='gray', c=optical_colors[band],label=band)


        if limit_dates!=[]:
            ax1.axvline(limit_dates[0], linestyle='--')
            ax1.axvline(limit_dates[1], linestyle='--')
        fig.canvas.mpl_connect('pick_event', lambda event: click_action(self.ra, self.dec, xmm_name, swift_name))

        order = np.argsort(time_track)
        hardness_track=np.array(hardness_track)[order]
        luminosity_track = np.array(luminosity_track)[order]
        hardness_err_track=[np.array(hardness_err_track[0])[order],np.array(hardness_err_track[1])[order]]
        luminosity_err_track=[np.array(luminosity_err_track[0])[order],np.array(luminosity_err_track[1])[order]]
        time_track = np.array(time_track)[order]
        catalogs_track = np.array(catalogs_track)[order]
        color_track = np.array([matplotlib.cm.get_cmap('inferno')((time-time_track[0])/(time_track[-1]-time_track[0])) for time in time_track])
        ax3.errorbar(hardness_track,luminosity_track,xerr=hardness_err_track,yerr=luminosity_err_track, alpha=0.2, linestyle="--")
        #ax3.scatter(hardness_track,luminosity_track, c=ax3.lines[-1].get_color())
        for cat in catalogs:
            ax3.scatter(hardness_track[catalogs_track==cat], luminosity_track[catalogs_track==cat], c=ax3.lines[-1].get_color(),marker=hr_track_markers[cat], s=50, label=cat, edgecolors="gray")

        #ax1.tick_params(axis='x', rotation=45)
        ax1.set_title("Long-term lightcurve (0.2-12 keV)")
        ax1.legend(loc="best")
        ax1.set_yscale('log')
        ax1.set_xlabel("Time (MJD)")
        ax1.set_ylabel(r"Flux ($erg.s^{-1}.cm^{-2}$)")
        margin = 0.1*(self.max_time-self.min_time)
        ax1.set_xlim(self.min_time-margin, self.max_time+margin)

        ax2.set_title("Detections X-ray spectra")
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.set_xlabel("Energy (keV)")
        ax2.set_ylabel(r"$F_{\nu}$ ($erg.s^{-1}.cm^{-2}.keV^{-1}$)")

        ax3.set_title("Hardness-Luminosity diagram")
        ax3.legend(loc="best")
        ax3.set_yscale('log')
        ax3.set_xlabel("Hardness")
        ax3.set_xlim((-1.1,1.1))
        ax3.set_ylabel(r"Flux ($erg.s^{-1}.cm^{-2}$)")

        if self.optical_sources!={}:
            ax4.set_title("Optical/UV lightcurves")
            ax4.set_yscale('log')
            ax4.set_xlabel("Time (MJD)")
            ax4.set_ylabel(r"Flux ($erg.s^{-1}.cm^{-2}$)")
            ax4.legend(loc="best")
            ax4.sharex(ax1)

            ax5.legend()
            ax5.set_title("SED")
            ax5.set_xlabel("Frequency (Hz)")
            ax5.set_ylabel(r'$\nu$ $F_{\nu}$ (Hz.erg.s$^{-1}.$cm$^{-2}$.Hz$^{-1}$)')
            ax5.set_xscale('log')
            ax5.set_yscale('log')

            ax6.set_title(rf"Time evolution of $\alpha$")
            ax6.set_xlabel("Time (MJD)")
            ax6.legend()
            ax6.sharex(ax1)


        #Finally, if there is a GLADE counterpart, we have distance so we convert all fluxes to luminosities, and
        #add an axis on the right of all concerned plots
        if self.glade_distance!=[] and self.flux_lum_conv_factor>0:
            second_axis_func = (lambda x:self.flux_lum_conv_factor*x, lambda x:x/self.flux_lum_conv_factor)
            for ax in [ax1, ax3]:
                secax = ax.secondary_yaxis('right', functions=second_axis_func)
                secax.set_ylabel(r'Luminosity ($erg.s^{-1})$')
            secax = ax2.secondary_yaxis('right', functions=second_axis_func)
            secax.set_ylabel(r'$L_{\nu}$ ($erg.s^{-1})$')
            if self.optical_sources != {}:
                secax = ax4.secondary_yaxis('right', functions=second_axis_func)
                secax.set_ylabel(r'Luminosity ($erg.s^{-1})$')
                secax = ax5.secondary_yaxis('right', functions=second_axis_func)
                secax.set_ylabel(r'$\nu$ $L_{\nu}$ (Hz.erg.s$^{-1}$.Hz$^{-1}$)')
        plt.draw()
        plt.tight_layout()
        plt.show()


def load_source(cat):
    """
    This loads the catalog data for a given X-ray Source catalog.
    :param cat: Name of the catalog, using the same naming convention as the "catalogs" Table.
    :return: Dictionary, with the name of the source as a key and the Source object as a value
    """
    print(f"Loading {cat}...")
    raw_data = fits.open(f"{path_to_catalogs}{cat}.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    sources_raw = sources_raw[np.argsort(sources_raw[src_names[cat]])]

    #This retrieves the indices in the catalog that correspond to a change of source
    indices_for_source = [i for i in range(1, len(sources_raw)) if (sources_raw[src_names[cat]][i] != sources_raw[src_names[cat]][i - 1])]

    if cat == "Swift":
        #For Swift, we need the Start and Stop times for the Stacked fluxes
        timestartobs = Time(sources_raw["StartTime_UTC"], format="iso").mjd
        timeendobs = Time(sources_raw["StopTime_UTC"], format="iso").mjd
        timestartobs = np.split(timestartobs, indices_for_source)
        timeendobs = np.split(timeendobs, indices_for_source)


    #We divide up the catalog in sub-samples corresponding to each source
    timesteps = np.split(np.array(sources_raw[time_names[cat]]), indices_for_source)
    if cat in ("XMM","Stacked"):
        xmm_offaxes = np.split(np.array(sources_raw["EP_OFFAX"]), indices_for_source)
    else:
        xmm_offaxes = [[] for elt in indices_for_source]
    if cat in ("XMM","Swift","Stacked"):
        obsids = np.split(np.array(sources_raw[obsid_names[cat]]), indices_for_source)
    else:
        obsids = [[] for elt in indices_for_source]
    if cat in ("XMM","Chandra"):
        short_term_var_flags = np.split(np.array(sources_raw[short_term_var_name[cat]]), indices_for_source)
    else:
        short_term_var_flags = [[] for elt in indices_for_source]
    names = np.split(np.array(sources_raw[src_names[cat]]), indices_for_source)

    band_fluxes = []
    band_flux_errors_neg=[]
    band_flux_errors_pos=[]

    fluxes = np.split(conv_factors[cat]*np.array(sources_raw[flux_names[cat]]), indices_for_source)
    flux_errors_neg = np.split(conv_factors[cat]*np.array(sources_raw[flux_error_names[cat][0]]), indices_for_source)
    flux_errors_pos = np.split(conv_factors[cat]*np.array(sources_raw[flux_error_names[cat][1]]), indices_for_source)
    flux_errors = [[flux_neg, flux_pos] for (flux_neg, flux_pos) in zip(flux_errors_neg, flux_errors_pos)]

    for band_flux_name, band_fluxerr_neg_name, band_fluxerr_pos_name in zip(band_flux_names[cat],
                                                                          band_fluxerr_names[cat][0],band_fluxerr_names[cat][1]):
        band_fluxes.append(np.array(sources_raw[band_flux_name]))
        band_flux_errors_neg.append(np.array(sources_raw[band_fluxerr_neg_name]))
        band_flux_errors_pos.append(np.array(sources_raw[band_fluxerr_pos_name]))

    #We transpose the band fluxes tables, so the first dimension corresponds to sources and not energy bands
    band_fluxes = np.transpose(np.array(band_fluxes))
    band_flux_errors_neg = np.transpose(np.array(band_flux_errors_neg))
    band_flux_errors_pos = np.transpose(np.array(band_flux_errors_pos))
    band_fluxes = np.split(band_fluxes, indices_for_source)
    band_flux_errors_neg = np.split(band_flux_errors_neg, indices_for_source)
    band_flux_errors_pos = np.split(band_flux_errors_pos, indices_for_source)
    band_fluxerr = [[band_flux_neg, band_flux_pos] for (band_flux_neg, band_flux_pos) in zip(band_flux_errors_neg, band_flux_errors_pos)]


    dic_sources = {}

    #This loops on all sources, to build the Source objects
    for (index, flux, flux_error, time, name, band_flux, band_fluxerr, obsid, xmm_offaxis, short_term_var) in zip(range(len(fluxes)),fluxes, flux_errors, timesteps, names, band_fluxes, band_fluxerr, obsids, xmm_offaxes, short_term_var_flags):
            swift_stacked_flux=[]
            swift_stacked_flux_err=[[],[]]
            swift_stacked_times=[[],[]]
            if cat == "Swift":
                tab_src_timestartobs = timestartobs[index]
                tab_src_timeendobs = timeendobs[index]

                #We select the stacked Swift detections first
                swift_stacked_flux=flux[obsid>1e10]
                swift_stacked_flux_err=[flux_error[0][obsid>1e10],flux_error[1][obsid>1e10]]
                swift_stacked_times=[tab_src_timestartobs[obsid>1e10], tab_src_timeendobs[obsid>1e10]]

                # We then treat the classical, non-stacked Swift detections
                flux = flux[obsid < 1e10]
                flux_error = [flux_error[0][obsid < 1e10], flux_error[1][obsid < 1e10]]
                time = time[np.where(obsid < 1e10)]
                band_flux = band_flux[obsid < 1e10]
                band_fluxerr = [band_fluxerr[0][obsid < 1e10], band_fluxerr[1][obsid < 1e10]]
                obsid = obsid[obsid < 1e10]
            source = Source(cat, name[0].strip(), flux, flux_error, time, band_flux, band_fluxerr, obsid, swift_stacked_flux,swift_stacked_flux_err,swift_stacked_times, xmm_offaxis, short_term_var)
            dic_sources[name[0].strip()] = source
    return dic_sources

def load_optical_source(cat):
    """
    This loads the catalog data for a given optical Source catalog.
    :param cat: Name of the catalog, here OM or UVOT
    :return: Dictionary, with the name of the source as a key and the OpticalSource object as a value
    """
    print(f"Loading {cat}...")
    raw_data = fits.open(f"{path_to_catalogs}{cat}.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    sources_raw = sources_raw[np.argsort(sources_raw[src_names[cat]])]

    indices_for_source = [i for i in range(1, len(sources_raw)) if (sources_raw[src_names[cat]][i] != sources_raw[src_names[cat]][i - 1])]

    #We divide up the catalog in sub-samples corresponding to each source
    timesteps = np.split(np.array(sources_raw[time_names[cat]]), indices_for_source)
    obsids = np.split(np.array(sources_raw[obsid_names[cat]]), indices_for_source)
    names = np.split(np.array(sources_raw[src_names[cat]]), indices_for_source)

    band_fluxes = []
    band_flux_errors=[]
    for band_flux_name, band_fluxerr_name, halfband_width in zip(band_flux_names[cat], band_fluxerr_names[cat], band_half_width[cat]):
        band_fluxes.append(np.array(sources_raw[band_flux_name])*2*halfband_width)
        band_flux_errors.append(np.array(sources_raw[band_fluxerr_name])*halfband_width)

    #We use the same framework for the optical filters as for the energy bands in X-ray, so we need to transpose so the
    #first dimension corresponds to sources and not filters
    band_fluxes = np.transpose(np.array(band_fluxes))
    band_flux_errors = np.transpose(np.array(band_flux_errors))
    band_fluxes = np.split(band_fluxes, indices_for_source)
    band_flux_errors = np.split(band_flux_errors, indices_for_source)

    dic_sources = {}

    #This loops on all sources, to build the Source objects
    pbar=tqdm(total=len(band_fluxes))
    for (index, time, name, band_flux, band_fluxerr, obsid) in zip(range(len(band_fluxes)), timesteps, names, band_fluxes, band_flux_errors, obsids):
        source = OpticalSource(cat, name[0].strip(), time, band_flux, band_fluxerr, obsid)
        dic_sources[name[0].strip()] = source
        pbar.update(1)
    pbar.close()
    return dic_sources

def load_GLADE():
    """
    Loads the GLADE galaxy catalog data
    :return: A dictionary, with the GLADE identifier as a key, and a Table containing Distance and Galaxy Stellar Mass as value
    """
    print(f"Loading GLADE+...")
    raw_data = fits.open(f"{path_to_catalogs}GLADE.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    glade_distances ={}
    for line in sources_raw:
        glade_distances[line["GLADE_IAUNAME"]]=[float(line["d_L"]),float(line["stellar_mass"])]
    return glade_distances

def compute_x_to_optical(dic_master_sources):
    """
    Computes the associations between X-ray and Optical detections, the alpha_OX. Is also used to estimate the correlation
    between the X-ray and the optical lightcurves, through the computation of a Pearson coefficient. Its errors are estimated
    through a Monte Carlo method.
    :param dic_master_sources: One of the produce of load_master_sources(), which will call this function
    :return: Nothing, the MasterSource objects are directly modified
    """
    print("Computing X-to-optical ratio evolution...")
    for ms in tqdm(list(dic_master_sources.values())):
        ms.alpha_band_x_Fnu = {}
        ms.alpha_band_x_Fnu_error = {}
        ms.alpha_band_x_timesteps = {}
        ms.pearson = {}
        ms.pearson_stderr = {}
        ms.pearson_percentileerr={}
        x_fluxdens={}
        x_fluxdenserr={}
        opt_fluxdens={}
        opt_fluxdenserr={}
        for band in ["UVW2", "UVM2", "UVW1", "U", "B", "V"]:
            ms.alpha_band_x_Fnu[band]=[]
            ms.alpha_band_x_Fnu_error[band] = []
            ms.alpha_band_x_timesteps[band] = []
            x_fluxdens[band] = []
            opt_fluxdens[band] = []
            x_fluxdenserr[band] = []
            opt_fluxdenserr[band] = []
        for (opt_cat, x_cat) in zip(["OM", "UVOT"], ['XMM', 'Swift']):
            if (opt_cat in ms.optical_sources.keys()) and (x_cat in ms.sources.keys()):
                x_obsid = ms.sources[x_cat].obsids
                x_flux_densities = {}
                x_flux_density_errors = {}
                x_times = {}
                for obs, flux, fluxerr_neg, fluxerr_pos, timestep in zip(x_obsid, ms.sources[x_cat].fluxes, ms.sources[x_cat].flux_errors[0],ms.sources[x_cat].flux_errors[1], ms.sources[x_cat].timesteps):
                    x_flux_densities[obs] = flux / xband_width
                    x_flux_density_errors[obs] = [fluxerr_neg / xband_width,fluxerr_pos/ xband_width]
                    x_times[obs] = timestep
                opt_obsid = ms.optical_sources[opt_cat].obsids
                opt_flux_densities = {}
                opt_flux_density_errors = {}

                for band_ind, band in enumerate(["UVW2", "UVM2", "UVW1", "U", "B", "V"]):
                    opt_flux_densities[band] = {}
                    opt_flux_density_errors[band] = {}
                    for obs, flux, fluxerr in zip(opt_obsid, ms.optical_sources[opt_cat].band_flux, ms.optical_sources[opt_cat].band_fluxerr):
                        opt_flux_densities[band][obs] = flux[band_ind] / (2*frequencies_half_width[cat][band_ind])
                        opt_flux_density_errors[band][obs] = fluxerr[band_ind] / (2*frequencies_half_width[cat][band_ind])
                    all_obsids = list(set(list(x_obsid) + list(opt_obsid)))
                    for obsid in all_obsids:
                        #We loop on all X-ray and Optical ObsID and only keep those that are common
                        if (obsid in x_flux_densities.keys()) and (obsid in opt_flux_densities[band].keys()):
                            if (not np.isnan(opt_flux_densities[band][obsid])) and (not np.isnan(x_flux_densities[obsid])):
                                x_fluxdens[band].append(x_flux_densities[obsid])
                                x_fluxdenserr[band].append(np.mean(x_flux_density_errors[obsid]))
                                opt_fluxdens[band].append(opt_flux_densities[band][obsid])
                                opt_fluxdenserr[band].append(opt_flux_density_errors[band][obsid])
                                alpha_band_x_Fnu = (np.log10(x_flux_densities[obsid]) - np.log10(opt_flux_densities[band][obsid])) / (np.log10(xband_average_frequency) - np.log10(frequencies[opt_cat][band_ind]))
                                ms.alpha_band_x_Fnu[band].append(alpha_band_x_Fnu)
                                ms.alpha_band_x_Fnu_error[band].append([np.log10((1+opt_flux_density_errors[band][obsid]/opt_flux_densities[band][obsid])/(1-x_flux_density_errors[obsid][0]/x_flux_densities[obsid]))/(np.log10(xband_average_frequency)-np.log10(frequencies[opt_cat][band_ind])),
                                                                        np.log10((1+x_flux_density_errors[obsid][1]/x_flux_densities[obsid])/(1-opt_flux_density_errors[band][obsid]/opt_flux_densities[band][obsid]))/(np.log10(xband_average_frequency)-np.log10(frequencies[opt_cat][band_ind]))])
                                ms.alpha_band_x_timesteps[band].append(x_times[obsid])

        ms.alpha_band_x_Fnu_var = {}
        ms.alpha_band_x_nuFnu_var = {}
        #200 might be a little small depending on the cases, but a lot more would slow down the computation by too much
        n_sample=200
        for band in ["UVW2", "UVM2", "UVW1", "U", "B", "V"]:
            if len(x_fluxdens[band]) > 1:
                sampling_lc1 = np.random.normal(x_fluxdens[band], x_fluxdenserr[band], (n_sample,len(x_fluxdens[band])))
                sampling_lc2 = np.random.normal(opt_fluxdens[band], opt_fluxdenserr[band], (n_sample,len(x_fluxdens[band])))
                tab_pearson_sampling = [pearsonr(sampled1, sampled2)[0] for sampled1, sampled2 in
                                        zip(sampling_lc1, sampling_lc2)]
                ms.pearson[band] = np.mean(tab_pearson_sampling)
                ms.pearson_stderr[band] = np.std(tab_pearson_sampling)
                ms.pearson_percentileerr[band] = np.diff(np.percentile(tab_pearson_sampling,(16,84)))
            else:
                ms.pearson[band] = np.nan
                ms.pearson_stderr[band] = np.nan
                ms.pearson_percentileerr[band] = np.nan
            ms.alpha_band_x_Fnu_var[band]=0
            ms.alpha_band_x_nuFnu_var[band]=0
            if len(ms.alpha_band_x_Fnu[band])>1:
                max_loweralpha = np.nanmax(ms.alpha_band_x_Fnu[band]-np.transpose(ms.alpha_band_x_Fnu_error[band])[0])
                min_upperalpha = np.nanmin(ms.alpha_band_x_Fnu[band]+np.transpose(ms.alpha_band_x_Fnu_error[band])[1])
                ms.alpha_band_x_Fnu_var[band] = max_loweralpha-min_upperalpha

def load_XMM_upperlimits(dic_master_sources):
    """
    Loads the pre-computed XMM-Newton upper limits, which use the RapidXMM framework and the MasterSource identifiers.
    :param dic_master_sources
    :return: Nothing; MasterSource objects are updated with the corresponding Upper Limits
    """
    raw_data = fits.open(f"{path_to_master_sources}Master_source_XMM_UpperLimits.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    xmm_pointed_ul = sources_raw[sources_raw["obstype"]=="pointed"]
    dates = Time(xmm_pointed_ul["MJD_Date"],format="mjd").mjd
    print("Loading XMM pointed upper limits on Master Sources...")
    pbar = tqdm(total=len(dates))
    for line, date in zip(xmm_pointed_ul, dates):
        ms = dic_master_sources[line["MS_ID"]]
        if "Stacked" in ms.sources.keys():
            #We don't use upper limits in the case of a Stacked clean detection, as we assume the filtering of the Stacked
            #catalog was of better quality than the automatic RapidXMM pipeline
            if int(line["obsid"]) not in ms.sources["Stacked"].obsids:
                #In the case of a simultaneous Stacked detection / RapidXMM upper limit, we take the Stacked detection
                ms.xmm_ul.append(line["ul_flux8_3sig"])
                ms.xmm_ul_dates.append(date)
                ms.xmm_ul_obsids.append(line["obsid"])
                ms.min_time = min(date, ms.min_time)
                ms.max_time = max(date, ms.max_time)
        else:
                ms.xmm_ul.append(line["ul_flux8_3sig"])
                ms.xmm_ul_dates.append(date)
                ms.xmm_ul_obsids.append(line["obsid"])
                ms.min_time = min(date, ms.min_time)
                ms.max_time = max(date, ms.max_time)
        pbar.update(1)
    pbar.close()

    print("Loading XMM slew upper limits on Master Sources...")
    xmm_slew_ul = sources_raw[sources_raw["obstype"] == 'slew   ']
    dates = Time(xmm_slew_ul["MJD_Date"], format="mjd").mjd
    pbar = tqdm(total=len(dates))
    for line, date in zip(xmm_slew_ul, dates):
        ms = dic_master_sources[line["MS_ID"]]
        ms.slew_ul.append(line["ul_flux8_3sig"])
        ms.slew_ul_dates.append(date)
        ms.slew_ul_obsids.append(line["obsid"])
        ms.min_time = min(date, ms.min_time)
        ms.max_time = max(date, ms.max_time)
        pbar.update(1)
    pbar.close()


    #Updating variabilities using upper limits
    tab_improvement_withxxm = []
    tab_improvement_withoutxxm = []
    to_plot = []
    improve=[]
    for ms in dic_master_sources.values():
        if len(ms.xmm_ul)>0:
            if ms.min_upper <1 and min(ms.xmm_ul) > 0:
                if "XMM" in ms.sources.keys() or "Stacked" in ms.sources.keys():
                    tab_improvement_withxxm.append(ms.min_upper / min(ms.xmm_ul))
                else:
                    tab_improvement_withoutxxm.append(ms.min_upper / min(ms.xmm_ul))
                    if (ms.min_upper / min(ms.xmm_ul)>10):
                        to_plot.append(ms)
                        improve.append(ms.min_upper / min(ms.xmm_ul))
        ms.var = ms.max_lower/ms.min_upper

def load_Chandra_upperlimits(dic_master_sources):
    """
    Loads the pre-computed Chandra upper limits, which use the MasterSource identifiers.
    :param dic_master_sources
    :return: Nothing; MasterSource objects are updated with the corresponding Upper Limits
    """
    print("Load Chandra upper limits...")
    raw_data = fits.open(f"{path_to_master_sources}Chandra_UL.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    for line in tqdm(sources_raw):
        ms = dic_master_sources[line["MS_ID"]]
        ms.chandra_ul.append(line["flux_aper_hilim_b"])
        ms.chandra_ul_dates.append(line["gti_mjd_obs"])


"""
def load_optical_only_sources(tab_optical_sources,glade_galaxies):
    print("Loading OM & UVOT sources with no X-ray counterparts...")
    raw_data = fits.open(f"{path_to_master_sources}Combined_Optical.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    tab_optical_only = []
    for line in tqdm(sources_raw):
        opt_source_for_ms=[]
        has_xray_counterpart = False
        for cat in optical_catalogs:
            if line[f"{cat}_IAUNAME"] in tab_optical_sources[cat].keys():
                opt_source = tab_optical_sources[cat][line[f"{cat}_IAUNAME"]]
                current_ra = line["RA_1"]
                current_dec = line["DEC_1"]
                if opt_source.master_source == []:
                    opt_source_for_ms.append(opt_source)
                else:
                    has_xray_counterpart=True
        if (len(opt_source_for_ms)>0) and (not has_xray_counterpart):
            ms = MasterSource(-1,[],current_ra,current_dec,2,opt_source_for_ms)
            if line["GLADE_IAUNAME"] != -2147483648:  # and multiwavelength_information:
                (ms.glade_distance, ms.glade_stellar_mass) = glade_galaxies[line["GLADE_IAUNAME"]]
                ms.flux_lum_conv_factor = 4 * np.pi * (ms.glade_distance * 3.086E+24) ** 2
            tab_optical_only.append(ms)
    return tab_optical_only
"""

def load_master_sources(multiwavelength_information=False):
    """
    This function will call all previous functions and build the MasterSource catalog, using all archival data.
    The user can decide whether they want to use the OM & UVOT data
    :param multiwavelength_information: Boolean True if you want to load OM & UVOT data, False otherwise.
    :return: a Dictionary, containing the MasterSources as values and their identifiers as key. Also returns a table
    containing the OpticalSource objects with no X-ray counterpart, but for now this is empty table.
    """
    tab_catalog_sources = {}
    for cat in catalogs[:-1]:
        tab_catalog_sources[cat] = load_source(cat)
    tab_optical_sources = {}
    if multiwavelength_information:
        for opt_cat in optical_catalogs:
            tab_optical_sources[opt_cat] = load_optical_source(opt_cat)
    glade_galaxies = load_GLADE()

    print(f"Loading Master Sources...")
    raw_data = fits.open(f"{path_to_master_sources}Master_source_enhanced_Simbad.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    dic_master_sources = {}
    tab_optical_only=[]
    known_TDE_MS_ids = Table(fits.open(f"{path_to_master_sources}Master_source_TDE.fits", memmap=True)[1].data)["MS_ID"]
    for line in tqdm(sources_raw):
        tab_sources_for_this_ms = []
        tab_optical_sources_for_this_ms = []
        for cat in catalogs[:-1]:
            if line[cat]!='':
                name=line[cat].strip()
                if name in tab_catalog_sources[cat].keys():
                    tab_sources_for_this_ms.append(tab_catalog_sources[cat][name])
        if multiwavelength_information:
            for cat in optical_catalogs:
                if line[cat]!='':
                    name=line[cat].strip()
                    if name in tab_optical_sources[cat].keys():
                        tab_optical_sources_for_this_ms.append(tab_optical_sources[cat][name])
        ms_id = line["MS_ID"]
        ms = MasterSource(ms_id, tab_sources_for_this_ms, line["MS_RA"], line["MS_DEC"], line["MS_POSERR"], tab_optical_sources_for_this_ms)
        ms.simbad_type = line["Simbad_type"]
        if ms.simbad_type=="Galaxy":
            ms.simbad_galaxy_type = line["main_type"]
        if line["GLADE_ID"] != -2147483648:# and multiwavelength_information:
            (ms.glade_distance,ms.glade_stellar_mass) = glade_galaxies[line["GLADE_ID"]]
            ms.flux_lum_conv_factor = 4*np.pi*(ms.glade_distance*3.086E+24)**2
        ms.is_known_tde = False
        if ms_id in known_TDE_MS_ids:
            ms.is_known_tde=True
        dic_master_sources[ms_id] = ms
    if multiwavelength_information:
        compute_x_to_optical(dic_master_sources)
        #tab_optical_only = load_optical_only_sources(tab_optical_sources,glade_galaxies)
    load_XMM_upperlimits(dic_master_sources)
    load_Chandra_upperlimits(dic_master_sources)
    for ms in dic_master_sources.values():
        ms.is_known_qpe = False
        if "XMM" in ms.sources.keys():
            if ms.sources["XMM"].name in ("4XMM J011908.6-341130", "4XMM J130200.1+274657"):
                ms.is_known_qpe = True
    print("Master sources loaded!")
    return dic_master_sources, tab_optical_only


