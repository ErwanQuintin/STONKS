"""
The Search for Transient Objects in New detections using Known Sources (STONKS) software is meant to allow automatic
comparison between the detections of a new XMM-Newton observation, and all the available archival X-ray data. It will send
out an alert and save the lightcurve and available data for all sources with a long-term variability over a factor of 5.
The updated framework is meant to be used with a single input position and flux state.

This requires the script to be run with the same file structure as on the GitHub, to ensure access to hand-made catalog
data.

The LoadSpecificMasterSources.py module loads the archival X-ray MasterSource object that matches the input position.

This software was developed as part of the XMM2ATHENA project. This project has received funding from the European
Union's Horizon 2020 research and innovation programme under grant agreement n°101004168, the XMM2ATHENA project.

Author: Erwan Quintin, erwan.quintin@irap.omp.eu
"""


import numpy as np
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import Longitude, Latitude, Angle
from astroquery.hips2fits import hips2fits
from tqdm import tqdm
from astropy.constants import c
import shlex
import subprocess
import os
import json
from dict_utils import DictUtils
from constants import PATHTO
import cmasher as cmr
from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})

style="bmh"
cmap_to_use= "cmr.torch"#"turbo"


#path_to_catalogs = "/home/erwan/Documents/PhD/LongTermVariability/LongTermVariability_Python/NewMatchMethod/CleanCatalogs/Catalogs/FullData/"
#path_to_master_sources = "/home/erwan/Documents/PhD/LongTermVariability/LongTermVariability_Python/NewMatchMethod/CleanCatalogs/MasterSource/"

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
                "RASS":1/0.41,
                "WGACAT":1/0.41}
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
# band_conv_factors_soft = {"XMM":0.35/0.35,
#                           "NewXMM":0.35/0.35,
#                           "Chandra":0.35/0.28,
#                           "Swift":0.35/0.34,
#                           "Slew":0.35/0.35,
#                           "eRosita":0.35/0.35,
#                           "RASS":0.35/0.35,
#                           "WGACAT":0.35/0.35,
#                           "Stacked":0.35/0.35}
# band_conv_factors_hard = {"XMM":0.65/0.65,
#                           "NewXMM":0.65/0.65,
#                           "Chandra":0.65/0.41,
#                           "Swift":0.65/0.56,
#                           "Slew":0.65/0.65,
#                           "eRosita":0.65/0.25,
#                           "RASS":np.nan,
#                           "WGACAT":np.nan,
#                           "Stacked":0.65/0.65}


#Values of hardness ratios at which the spectral assumption might break down, leading to a warning in the alert
dic_soft_threshold={'XMM': -0.415,
                    'NewXMM': -0.415,
                    'Slew': -0.415,
                    'Stacked': -0.415,
                    'Chandra': -0.33,
                    'Swift': -0.4,
                    'eRosita': -0.62}#For Gamma=2.5
dic_hard_threshold={'XMM': 0.88,
                    'NewXMM': 0.88,
                    'Slew': 0.88,
                    'Stacked': 0.88,
                    'Chandra': 0.74,
                    'Swift': 0.84,
                    'eRosita': 0.45}#For Gamma=0.5

hr_track_markers = {"XMM":"o","NewXMM":'o',"Chandra":"v","Swift":"s","Slew":"P","eRosita":"^","RASS":"d","WGACAT":"d","Stacked":"*"}

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

        self.soft_dets = [np.sum(det[:hr_bandlimit_index[catalog]]) for det in self.band_flux]
        self.soft_errors = [[np.sqrt(np.sum(np.array(err_neg[:hr_bandlimit_index[catalog]])**2)) for err_neg in self.band_fluxerr[0]],
                            [np.sqrt(np.sum(np.array(err_pos[:hr_bandlimit_index[catalog]])**2)) for err_pos in self.band_fluxerr[1]]]
        if catalog!= "RASS" and catalog!="WGACAT":
            self.hard_dets = [np.sum(det[hr_bandlimit_index[catalog]:]) for det in self.band_flux]
            self.hard_errors = [
                [np.sqrt(np.sum(np.array(err_neg[hr_bandlimit_index[catalog]:])**2)) for err_neg in
                 self.band_fluxerr[0]],
                [np.sqrt(np.sum(np.array(err_pos[hr_bandlimit_index[catalog]:])**2)) for err_pos in
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
        self.has_spectral_warning = False
        notnanhardness = np.array(self.hardness)[~np.isnan(self.hardness)]
        if len(notnanhardness)>0:
            self.has_spectral_warning = False in ((np.array(self.hardness)<dic_hard_threshold[catalog])
                                                  &(np.array(self.hardness)>dic_soft_threshold[catalog]))
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
        if len(swift_stacked_flux) > 0 :
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
    def __init__(self, session, id, tab_sources, ra, dec, poserr, tab_optical_sources):
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
        self.session = session
        self.id = id
        self.sources = {}
        self.sources_fluxes = []
        self.sources_error_bars = []
        self.sources_timesteps = []
        self.sources_var = []
        self.tab_hr = []
        self.tab_hr_err = [[],[]]
        self.never_on_axis_xmm = False
        self.has_short_term_var = False
        self.min_time=60000
        self.max_time=0
        self.has_spectral_warning=False
        for source in tab_sources:
            if ("XMM" in self.sources.keys()) and (source.catalog == "Stacked"):
                #We remove the Stacked detection that correspond to a clean XMM detection
                xmm_obsid = self.sources["XMM"].obsids
                stacked_obsid = source.obsids
                new_det_ind = [i for i in range(len(stacked_obsid)) if stacked_obsid[i] not in xmm_obsid]
                source.fluxes = source.fluxes[new_det_ind]
                source.flux_errors[0] = np.array(source.flux_errors[0])[new_det_ind]
                source.flux_errors[1] = np.array(source.flux_errors[1])[new_det_ind]
                source.timesteps = np.array(source.timesteps)[new_det_ind]
                source.obsids = np.array(source.obsids)[new_det_ind]
                source.hardness = np.array(source.hardness)[new_det_ind]
                source.hardness_err[0] = np.array(source.hardness_err[0])[new_det_ind]
                source.hardness_err[1] = np.array(source.hardness_err[1])[new_det_ind]

                source.band_flux = source.band_flux[new_det_ind]
                source.band_fluxerr[0] = np.array(source.band_fluxerr[0])[new_det_ind]
                source.band_fluxerr[1] = np.array(source.band_fluxerr[1])[new_det_ind]
            source.master_source = self
            self.sources[source.catalog]=source
            for (flux, fluxerr_neg, fluxerr_pos, timestep) in zip(source.fluxes, source.flux_errors[0], source.flux_errors[1], source.timesteps):
                self.sources_fluxes.append(flux)
                self.sources_error_bars.append(max(fluxerr_neg, fluxerr_pos))
                self.sources_var.append(source.var)
                self.sources_timesteps.append(timestep)
            self.tab_hr += list(source.hardness)
            self.tab_hr_err[0] += list(source.hardness_err[0])
            self.tab_hr_err[1] += list(source.hardness_err[1])
            for (flux, fluxerr_neg, fluxerr_pos, start, stop) in zip(source.swift_stacked_flux, source.swift_stacked_flux_err[0], source.swift_stacked_flux_err[1], source.swift_stacked_times[0], source.swift_stacked_times[1]):
                self.sources_fluxes.append(flux)
                self.sources_error_bars.append(max(fluxerr_neg, fluxerr_pos))
                self.min_time = min(start, self.min_time)
                self.max_time = max(stop, self.max_time)
                self.sources_timesteps.append((start+stop)/2)
            if len(source.xmm_offaxis) != 0:
                if np.nanmin(source.xmm_offaxis)>1:
                    self.never_on_axis_xmm = True
            if len(source.timesteps) != 0:
                self.min_time = min(min(source.timesteps), self.min_time)
                self.max_time = max(max(source.timesteps), self.max_time)
            for var_flag in source.short_term_var:
                if var_flag>0:
                    self.has_short_term_var=True
            if source.has_spectral_warning:
                self.has_spectral_warning=True
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
        self.simbad_name=''


    def save_lightcurve(self,dict_new_det_info, flag_alert, image_data, image_wcs):
        """
        Produces a multi-panel plot with most of the useful multi-instrument information about the source, saving it in
         a PDF file. From left to right and top to bottom:
        1. Long term multi-instrument X-ray lightcurves, extrapolated to the 0.1-12 keV
        2. Multi-instrument X-ray spectra, used to assess if the apparent flux change are due to a wrong extrapolation
        (if the spectra are comparable in common bands, but extrapolated fluxes are not)
        3. DSS Image
        :return: Nothing
        """
        cmap = cmr.take_cmap_colors('cmr.ocean', N=8, cmap_range=(0., 0.9))[::-1]
        colors = {}
        order = [7, 0, 4, 1, 3, 2, 5, 6]
        for cat, ind in zip(catalogs, order):
            colors[cat] = cmap[ind]
        colors["NewXMM"] = "darkred"

        fig, [[ax1, ax2], [ax3,ax4]] = plt.subplots(2,2, figsize=(10,10))


        if len(self.xmm_ul)!= 0:
            ax1.errorbar(Time(self.xmm_ul_dates,format='mjd').decimalyear, self.xmm_ul, yerr=0.2 * np.array(self.xmm_ul), uplims=True, fmt='none', c=colors["XMM"], label="XMM non-det.")
        if len(self.slew_ul)!= 0:
            ax1.errorbar(Time(self.slew_ul_dates,format='mjd').decimalyear, self.slew_ul, yerr=0.2 * np.array(self.slew_ul), uplims=True, fmt='none', c=colors["Slew"], label="Slew non-det.")
        if len(self.chandra_ul)!= 0:
            ax1.errorbar(Time(self.chandra_ul_dates,format='mjd').decimalyear, self.chandra_ul, yerr=0.2 * np.array(self.chandra_ul), uplims=True, fmt='none', c=colors["Chandra"])

        for cat in catalogs:
            if cat in self.sources.keys():
                source = self.sources[cat]
                #Plot the X-ray full band lightcurves
                ax1.errorbar(Time(np.array(source.timesteps),format='mjd').decimalyear, np.array(source.fluxes),
                            yerr=np.array(source.flux_errors), fmt="o",c=colors[cat],
                            ms=10,label=source.name, markeredgecolor='gray')
                if cat == "Swift":
                    times=[(stop+start)/2 for (start,stop) in zip(source.swift_stacked_times[0],source.swift_stacked_times[1])]
                    timerange=[(stop-start)/(2*365) for (start,stop) in zip(source.swift_stacked_times[0],source.swift_stacked_times[1])]
                    ax1.errorbar(Time(times,format='mjd').decimalyear, source.swift_stacked_flux,
                                yerr=source.swift_stacked_flux_err, xerr=timerange,
                                 fmt="o", markeredgecolor='gray', c=colors[cat])
                tab_width = 2*np.array(band_half_width[cat])
                for det in range(len(source.band_flux)):
                    # Plot the X-ray spectra of each detection, using the fluxes in various bands
                    good_indices=np.where(np.array(source.band_flux[det])>0)
                    ax2.plot(np.array(band_center[cat])[good_indices], np.array(source.band_flux[det])[good_indices] / tab_width[good_indices], c=colors[cat], lw=3)
                    # make sure the errors to be positive
                    yerr0 = [abs(number) for number in source.band_fluxerr[0][det]]
                    yerr1 = [abs(number) for number in source.band_fluxerr[1][det]]
                    ax2.errorbar(band_center[cat], source.band_flux[det] / tab_width,
                                 yerr=[yerr0 / tab_width,
                                       yerr1 / tab_width],
                                 xerr=band_half_width[cat],
                                 fmt="", c=colors[cat], alpha=0.2)
                    ax2.scatter(band_center[cat], source.band_flux[det] / tab_width, facecolor=colors[cat], marker="o",
                                edgecolor='gray', zorder=1)
                ax2.errorbar([],[],[],[],fmt='o', markeredgecolor='gray', c=colors[cat],label=cat)

        #Plot the X-ray image
        skyfov = 4 * u.arcmin
        skyposition = SkyCoord(self.ra * u.deg, self.dec * u.deg)
        pixposition = wcs.utils.skycoord_to_pixel(skyposition, image_wcs)
        pixelsize = wcs.utils.proj_plane_pixel_scales(image_wcs)[0]
        pixelfov = (skyfov.to(u.deg) / pixelsize).value

        imagecropped = image_data[int(pixposition[1] - (pixelfov / 2)):int(pixposition[1] + (pixelfov / 2)),
                                  int(pixposition[0] - (pixelfov / 2)):int(pixposition[0] + (pixelfov / 2))]
        ax3.imshow(np.where(imagecropped > 0, imagecropped, 1), norm=LogNorm(),cmap='cmr.ocean',origin='lower')
        ax3.scatter([pixelfov / 2], [pixelfov / 2], marker='+', c='r')
        ax3.axis("off")

        x0, y0 = 0, 0
        positions_x = [x0 + (0.25 * pixelfov), x0 + (0.5 * pixelfov)]
        positions_y = [y0 + (0.15 * pixelfov), y0 + (0.15 * pixelfov)]
        text_position_x = x0 + (0.38 * pixelfov)
        text_position_y = y0 + (0.07 * pixelfov)
        ax3.plot(positions_x, positions_y, c="w",lw=3)
        ax3.scatter(positions_x, positions_y, c="w", marker="o", s=20)
        scaletext = int(skyfov.value/4)
        ax3.text(text_position_x, text_position_y, f"{scaletext}'", c="w", fontsize=20,
                 horizontalalignment='center')

        #ax1.tick_params(axis='x', rotation=45)
        ax1.set_title("Long-term lightcurve (0.2-12 keV)")
        ax1.legend()
        ax1.set_yscale('log')
        ax1.set_xlabel("Time")
        ax1.set_ylabel(r"Flux ($erg.s^{-1}.cm^{-2}$)")
        ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
        #margin = 0.1*(self.max_time-self.min_time)
        #ax1.set_xlim(self.min_time-margin, self.max_time+margin)

        ax2.set_title("Detections X-ray spectra")
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.set_xlabel("Energy (keV)")
        ax2.set_ylabel(r"$F_{\nu}$ ($erg.s^{-1}.cm^{-2}.keV^{-1}$)")
        ax2.legend()


        #Finally, if there is a GLADE counterpart, we have distance so we convert all fluxes to luminosities, and
        #add an axis on the right of all concerned plots
        if self.glade_distance!=[] and self.flux_lum_conv_factor>0:
            second_axis_func = (lambda x:self.flux_lum_conv_factor*x, lambda x:x/self.flux_lum_conv_factor)
            secax = ax1.secondary_yaxis('right', functions=second_axis_func)
            secax.set_ylabel(r'Luminosity ($erg.s^{-1})$')
            secax = ax2.secondary_yaxis('right', functions=second_axis_func)
            secax.set_ylabel(r'$L_{\nu}$ ($erg.s^{-1})$')

        #Adding the metadata on the fourth panel
        obsid= DictUtils.get_value_by_key(dict_new_det_info, "ObsID")
        scr_num =  DictUtils.get_value_by_key(dict_new_det_info, "SRCNUM")

        ax4.set_xlim((0,1))
        ax4.set_ylim((0,1))
        r = fig.canvas.get_renderer()


        t1 = ax4.text(0., 0.92,  r"\textbf{Angular separation from target}", fontweight='bold')
        bb1 = t1.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        t2 = ax4.text(1., 0.92, str(DictUtils.get_value_by_key(dict_new_det_info, 'Angular separation from target')), horizontalalignment='right')
        bb2 = t2.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        ax4.plot([0. + bb1.width, 1. - bb2.width], [0.92 + 0.01, 0.92 + 0.01], ls=':', c='k')

        t1 = ax4.text(0., 0.87,  r"\textbf{Can be made public}", fontweight='bold')
        bb1 = t1.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        t2 = ax4.text(1., 0.87, str(DictUtils.get_value_by_key(dict_new_det_info, 'Publishable')), horizontalalignment='right')
        bb2 = t2.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        ax4.plot([0. + bb1.width, 1. - bb2.width], [0.87 + 0.01, 0.87 + 0.01], ls=':', c='k')


        for posy, key in zip([0.8, 0.75, 0.7],['ObsID','Date Obs','Exposure Time']):
            t1 = ax4.text(0., posy, r'\textbf{'+key+'}', fontweight='bold')
            bb1 = t1.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
            text_value = str(DictUtils.get_value_by_key(dict_new_det_info, key))
            # escape characters not supported
            text_value = text_value.replace("&", "\\&").replace("_", "\\_")
            t2 = ax4.text(1., posy,  text_value, horizontalalignment='right')
            bb2 = t2.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
            ax4.plot([0. + bb1.width, 1. - bb2.width], [posy+0.01, posy+0.01], ls=':', c='k')

        for posy, key in zip(np.arange(0.45,0.65, 0.05)[::-1],['SRCNUM', 'Source RA', 'Source Dec', 'Position Error']):
            t1 = ax4.text(0., posy, r'\textbf{'+key+'}', fontweight='bold')
            bb1 = t1.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
            t2 = ax4.text(1., posy, str(DictUtils.get_value_by_key(dict_new_det_info, key)), horizontalalignment='right')
            bb2 = t2.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
            ax4.plot([0. + bb1.width, 1. - bb2.width], [posy+0.01, posy+0.01], ls=':', c='k')
        t1=ax4.text(0., 0.4, r'\textbf{Instruments DetML}', fontweight='bold')
        bb1 = t1.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        t2=ax4.text(1., 0.35, str(DictUtils.get_value_by_key(dict_new_det_info, 'Instruments DetML')), horizontalalignment='right')
        bb2 = t2.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        ax4.plot([0. + bb1.width, 1.], [0.4 + 0.01, 0.4 + 0.01], ls=':', c='k')

        t1 = ax4.text(0., 0.25, r"\textbf{Type of Alert}", fontweight='bold')
        bb1 = t1.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        t2 = ax4.text(1., 0.25, flag_alert[0], horizontalalignment='right')
        bb2 = t2.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        ax4.plot([0. + bb1.width, 1. - bb2.width], [0.25 + 0.01, 0.25 + 0.01], ls=':', c='k')

        t1 = ax4.text(0., 0.2, r"\textbf{Long-term Variability}", fontweight='bold')
        bb1 = t1.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        t2 = ax4.text(1., 0.2, f'{np.round(self.var_ratio,1)}', horizontalalignment='right')
        bb2 = t2.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        ax4.plot([0. + bb1.width, 1. - bb2.width], [0.2 + 0.01, 0.2 + 0.01], ls=':', c='k')

        t1 = ax4.text(0., 0.15,  r"\textbf{Short-term Variability}", fontweight='bold')
        bb1 = t1.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        t2 = ax4.text(1., 0.15, f'{self.has_short_term_var}', horizontalalignment='right')
        bb2 = t2.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
        ax4.plot([0. + bb1.width, 1. - bb2.width], [0.15 + 0.01, 0.15 + 0.01], ls=':', c='k')

        if self.simbad_type.strip() !='':
            t1 = ax4.text(0., 0.1,  r"\textbf{Simbad}", fontweight='bold')
            bb1 = t1.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
            t2 = ax4.text(1., 0.1, f'{self.simbad_type.strip()} ({self.simbad_name.strip()})', horizontalalignment='right')
            bb2 = t2.get_window_extent(renderer=r).transformed(ax4.transData.inverted())
            ax4.plot([0. + bb1.width, 1. - bb2.width], [0.1 + 0.01, 0.1 + 0.01], ls=':', c='k')
        else:
            ax4.text(0., 0.1,  r"\textbf{Not in Simbad}", fontweight='bold')

        if self.has_spectral_warning:
            t1 = ax4.text(0., 0.05, r"\textbf{!!! Warning: extreme spectrum might impact variability !!!}", fontweight='bold')

        ax4.axis('off')
        plt.draw()
        plt.tight_layout()
        if obsid not in os.listdir(os.path.join(PATHTO.master_sources,'AlertsLightcurves')):
            os.mkdir(os.path.join(PATHTO.master_sources,'AlertsLightcurves',obsid))
        if scr_num.isdigit() is True:
            str_num = "{:03x}".format(int(scr_num))
        else:
            str_num = "FFF"
        filename = f"P{obsid}CAX000VALERT0{str_num}.PDF"
        lc_path = os.path.join(self.session.path,
                               filename)
        print(f"Saving Lightcurve {lc_path}")
        plt.savefig(lc_path)

    def save_json_alert(self, dict_new_det_info, flag_alert, param_holder):
        """Save the alert properties as a JSON file, to be used for a database"""
        obsid= DictUtils.get_value_by_key(dict_new_det_info, "ObsID")
        scr_num =  DictUtils.get_value_by_key(dict_new_det_info, "SRCNUM")
        if obsid not in os.listdir(os.path.join(PATHTO.master_sources,'AlertsLightcurves')):
            os.mkdir(os.path.join(PATHTO.master_sources,'AlertsLightcurves',obsid))
        if scr_num.isdigit() is True:
            str_num = "{:03x}".format(int(scr_num))
        else:
            str_num = "FFF"

        # Edit the detection infos from AlertPDF-oriented to JSON-oriented
        dict_new_det_info.pop("Off-axis Angles", None)
        dict_new_det_info.pop("Instruments DetML", None)
        dict_new_det_info.pop("Target angular separation", None)
        dict_new_det_info['EP_DETML']=f'{param_holder.ep_detml:.1f}'
        dict_new_det_info['PN_DETML']=f'{param_holder.pn_detml:.1f}'
        dict_new_det_info['M1_DETML']=f'{param_holder.m1_detml:.1f}'
        dict_new_det_info['M2_DETML']=f'{param_holder.m2_detml:.1f}'
        #dict_new_det_info['PN_OFFAX']=f'{param_holder.pn_offax:.1f}'
        #dict_new_det_info['M1_OFFAX']=f'{param_holder.m1_offax:.1f}'
        #dict_new_det_info['M2_OFFAX']=f'{param_holder.m2_offax:.1f}'
        dict_new_det_info['TARGET_ANG_SEP']=f'{param_holder.off_target_angle:.1f}'

        dict_new_det_info['VarAmplitude']=np.round(self.var_ratio,1)
        dict_new_det_info['LastFlux']=f'{param_holder.flux:.1e}'
        dict_new_det_info['LastFluxErr'] = f'{param_holder.flux_err:.1e}'
        last_soft, last_hard = np.nansum(param_holder.band_flux[0][:3]), np.nansum(param_holder.band_flux[0][3:])
        dict_new_det_info['LastHR']=f'{(last_hard-last_soft)/(last_hard+last_soft):.2f}'
        dict_new_det_info['ArchivalShortTermVar']=f'{self.has_short_term_var}'
        if self.simbad_type.strip() !='':
            dict_new_det_info['Simbad']= f'{self.simbad_type.strip()} ({self.simbad_name.strip()})'
        else:
            dict_new_det_info['Simbad']= 'Not in Simbad'

        dict_new_det_info['SpectralWarning']=f'{self.has_spectral_warning}'

        filename = f"P{obsid}CAX000VALERT0{str_num}.json"
        json_path = os.path.join(self.session.path, filename)
        dict_new_det_info['AlertType']=flag_alert[0]
        with open(json_path, "w") as outfile:
            json.dump(dict_new_det_info, outfile, indent=1)


def load_source_on_position(session, cat, ra_target, dec_target):
    """
    This loads the catalog data for the sources around a given position. Starts by constraining the catalog to the
    position, then loads the corresponding Source objects.

    :return: Dictionary, with the name of the source as a key and the Source object as a value
    """
    #print(f"Loading {cat}...")
    cmd = f"{PATHTO.stilts_cmd} tpipe {os.path.join(PATHTO.catalogs,str(cat)+'.fits')} cmd='select \"skyDistanceDegrees(RA,DEC,{ra_target},{dec_target})*60<20\"' \
    out={os.path.join(session.paths, str(cat)+'_MatchOnSource.fits')}"
    cmd = shlex.split(cmd)
    subprocess.run(cmd)
    raw_data = fits.open(os.path.join(session.path, str(cat)+'_MatchOnSource.fits'), memmap=True)
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
    if not (len(fluxes)==1 and list(fluxes[0])==[]):
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
    #cmd_cleanup=f"rm {PATHTO.catalogs}{cat}_MatchOnSource.fits"
    #cmd_cleanup = shlex.split(cmd_cleanup)
    #subprocess.run(cmd_cleanup)
    return dic_sources

def load_source_on_name(session, cat, given_name, ra_target, dec_target):
    """
    This loads the catalog data for the sources around a given position. Starts by constraining the catalog to the
    position, then loads the corresponding Source objects.
    :param cat: Name of the catalog, using the same naming convention as the "catalogs" Table.
    :return: Dictionary, with the name of the source as a key and the Source object as a value
    """
    #print(f"Loading {cat}...")
    cmd = f"{PATHTO.stilts_cmd} tpipe {os.path.join(PATHTO.catalogs,str(cat)+'.fits')} cmd='select \"skyDistanceDegrees(RA,DEC,{ra_target},{dec_target})*60<1\"' \
    out={os.path.join(session.path,str(cat)+'_MatchOnSource.fits')}"
    cmd = shlex.split(cmd)
    subprocess.run(cmd)
    raw_data = fits.open(os.path.join(session.path,str(cat)+'_MatchOnSource.fits'), memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    sources_raw = sources_raw[sources_raw[src_names[cat]]==given_name]

    source = []

    if len(sources_raw)>0:
        if cat == "Swift":
            #For Swift, we need the Start and Stop times for the Stacked fluxes
            timestartobs = Time(sources_raw["StartTime_UTC"], format="iso").mjd
            timeendobs = Time(sources_raw["StopTime_UTC"], format="iso").mjd

        #We divide up the catalog in sub-samples corresponding to each source
        timesteps = np.array(sources_raw[time_names[cat]])
        if cat in ("XMM","Stacked"):
            xmm_offaxes = np.array(sources_raw["EP_OFFAX"])
        else:
            xmm_offaxes = []
        if cat in ("XMM","Swift","Stacked"):
            obsids = np.array(sources_raw[obsid_names[cat]])
        else:
            obsids = []
        if cat in ("XMM","Chandra"):
            short_term_var_flags = np.array(sources_raw[short_term_var_name[cat]])
        else:
            short_term_var_flags = []

        band_fluxes = []
        band_flux_errors_neg=[]
        band_flux_errors_pos=[]

        fluxes = conv_factors[cat]*np.array(sources_raw[flux_names[cat]])
        flux_errors_neg = conv_factors[cat]*np.array(sources_raw[flux_error_names[cat][0]])
        flux_errors_pos = conv_factors[cat]*np.array(sources_raw[flux_error_names[cat][1]])
        flux_errors = [flux_errors_neg, flux_errors_pos]

        for band_flux_name, band_fluxerr_neg_name, band_fluxerr_pos_name in zip(band_flux_names[cat],
                                                                              band_fluxerr_names[cat][0],band_fluxerr_names[cat][1]):
            band_fluxes.append(np.array(sources_raw[band_flux_name]))
            band_flux_errors_neg.append(np.array(sources_raw[band_fluxerr_neg_name]))
            band_flux_errors_pos.append(np.array(sources_raw[band_fluxerr_pos_name]))

        #We transpose the band fluxes tables, so the first dimension corresponds to sources and not energy bands
        band_fluxes = np.transpose(np.array(band_fluxes))
        band_flux_errors_neg = np.transpose(np.array(band_flux_errors_neg))
        band_flux_errors_pos = np.transpose(np.array(band_flux_errors_pos))
        band_fluxerr = [band_flux_errors_neg,band_flux_errors_pos]

        #This loops on all sources, to build the Source objects
        if True:#not (len(fluxes)==1 and fluxes==[]):
            swift_stacked_flux=[]
            swift_stacked_flux_err=[[],[]]
            swift_stacked_times=[[],[]]
            if cat == "Swift":
                #We select the stacked Swift detections first
                swift_stacked_flux=fluxes[obsids>1e10]
                swift_stacked_flux_err=[flux_errors[0][obsids>1e10],flux_errors[1][obsids>1e10]]
                swift_stacked_times=[timestartobs[obsids>1e10], timeendobs[obsids>1e10]]

                # We then treat the classical, non-stacked Swift detections
                fluxes = fluxes[obsids < 1e10]
                flux_errors = [flux_errors[0][obsids < 1e10], flux_errors[1][obsids < 1e10]]
                timesteps = timesteps[np.where(obsids < 1e10)]
                band_fluxes = band_fluxes[obsids < 1e10]
                band_fluxerr = [band_fluxerr[0][obsids < 1e10], band_fluxerr[1][obsids < 1e10]]
                obsids = obsids[obsids < 1e10]
            source = Source(cat, given_name, fluxes, flux_errors, timesteps, band_fluxes, band_fluxerr, obsids,
                                swift_stacked_flux,swift_stacked_flux_err,swift_stacked_times, xmm_offaxes, short_term_var_flags)
    #cmd_cleanup=f"rm {PATHTO.catalogs}{cat}_MatchOnSource.fits"
    #cmd_cleanup = shlex.split(cmd_cleanup)
    #subprocess.run(cmd_cleanup)
    return source

def load_GLADE():
    """
    Loads the GLADE galaxy catalog data
    :return: A dictionary, with the GLADE identifier as a key, and a Table containing Distance and Galaxy Stellar Mass as value
    """
    print(f"Loading GLADE+...")
    raw_data = fits.open(os.path.join(PATHTO.catalogs,"GLADE.fits"), memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    glade_distances ={}
    for line in sources_raw:
        glade_distances[line["GLADE_IAUNAME"]]=[float(line["d_L"]),float(line["stellar_mass"])]
    return glade_distances




def load_specific_XMM_upperlimits(dic_master_sources, ms_id, obsid):
    """
    Loads the pre-computed XMM-Newton upper limits, which use the RapidXMM framework and the MasterSource identifiers.
    :param dic_master_sources
    :return: Nothing; MasterSource objects are updated with the corresponding Upper Limits
    """
    raw_data = fits.open(os.path.join(PATHTO.precomputed_obsids,'UpperLimits_'+str(obsid)+'.fits'), memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    sources_raw = sources_raw[sources_raw["MS_ID"]==ms_id]

    xmm_pointed_ul = sources_raw[sources_raw["obstype"] == "pointed"]
    dates = Time(xmm_pointed_ul["MJD_Date"],format="mjd").mjd
    ms = dic_master_sources[ms_id]
    for line, date in zip(xmm_pointed_ul, dates):
        if "Stacked" in ms.sources.keys():
            #We don't use upper limits in the case of a Stacked clean detection, as we assume the filtering of the Stacked
            #catalog was of better quality than the automatic RapidXMM pipeline
            if int(line["obsid"]) not in ms.sources["Stacked"].obsids:
                #In the case of a simultaneous Stacked detection / RapidXMM upper limit, we take the Stacked detection
                ms.xmm_ul.append(line["ul_flux8_3sig"])
                ms.xmm_ul_dates.append(date)
                ms.xmm_ul_obsids.append(line["obsid"])
        else:
                ms.xmm_ul.append(line["ul_flux8_3sig"])
                ms.xmm_ul_dates.append(date)
                ms.xmm_ul_obsids.append(line["obsid"])

    xmm_slew_ul = sources_raw[sources_raw["obstype"] == 'slew   ']
    dates = Time(xmm_slew_ul["MJD_Date"], format="mjd").mjd
    for line, date in zip(xmm_slew_ul, dates):
        ms.slew_ul.append(line["ul_flux8_3sig"])
        ms.slew_ul_dates.append(date)
        ms.slew_ul_obsids.append(line["obsid"])

def load_Chandra_upperlimits(dic_master_sources):
    """
    Loads the pre-computed Chandra upper limits, which use the MasterSource identifiers.
    :param dic_master_sources
    :return: Nothing; MasterSource objects are updated with the corresponding Upper Limits
    """
    print("Load Chandra upper limits...")
    raw_data = fits.open(os.path.join(PATHTO.master_sources,"Chandra_UL.fits"), memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    for line in tqdm(sources_raw):
        ms = dic_master_sources[line["MS_ID"]]
        ms.chandra_ul.append(line["flux_aper_hilim_b"])
        ms.chandra_ul_dates.append(line["gti_mjd_obs"])

def load_specific_master_sources(session, ms_id,obsid, ra_target, dec_target):
    """
    This function will call all previous functions and build the MasterSource, using all relevant archival data.
    For now, we don't load the OM & UVOT data
    :param ms_id, and obsid, which are the ID of the MasterSource and the Observation currently considered. The RA
    and Dec of the target are also given. These parameters are meant to allow pre-computations of sub-regions of the
    archival catalog, to increase speed.
    :return: a Dictionary, containing the MasterSource as value and the identifier as key.
    """
    raw_data = fits.open(os.path.join(PATHTO.master_sources,'PreComputedObsidMatches',str(obsid)+'.fits'), memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    line = sources_raw[sources_raw["MS_ID"]==ms_id]
    tab_sources_for_this_ms = []
    dic_master_sources={}
    for cat in catalogs[:-1]:
        if line[cat]!='':
            name=line[cat][0].strip()
            source = load_source_on_name(session, cat, name,  ra_target, dec_target)
            if source != []:
                tab_sources_for_this_ms.append(source)
        ms_id = line["MS_ID"][0]
        ms = MasterSource(session, ms_id, tab_sources_for_this_ms, line["MS_RA"], line["MS_DEC"], line["MS_POSERR"],[])
        ms.simbad_type = line["Simbad_type"]
        if ms.simbad_type=="Galaxy":
            ms.simbad_galaxy_type = line["main_type"]
        #if line["GLADE_ID"] != -2147483648:# and multiwavelength_information:
        #    (ms.glade_distance,ms.glade_stellar_mass) = glade_galaxies[line["GLADE_ID"]]
        #    ms.flux_lum_conv_factor = 4*np.pi*(ms.glade_distance*3.086E+24)**2
        dic_master_sources[ms_id] = ms
    load_specific_XMM_upperlimits(dic_master_sources,ms_id,obsid)
    #load_Chandra_upperlimits(dic_master_sources)
    #print("Master sources loaded!")
    #cmd_cleanup=f"rm {path_to_master_sources}Master_source_MatchOnSource.fits"
    #cmd_cleanup = shlex.split(cmd_cleanup)
    #subprocess.run(cmd_cleanup)
    return dic_master_sources

