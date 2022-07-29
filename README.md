The Search for Transient Objects in New detections using Known Sources (STONKS) software is meant to allow automatic
comparison between the PPS file of a new XMM-Newton observation, and all the available archival X-ray data. It will send
out an alert and plot the lightcurve and available data for all sources with a long-term variability over a factor of 5.

This repository contains 3 Python scripts:
1. **LoadMasterSources.py**, which is used to load the archival X-ray data in the form of MasterSource objects;
2. **STONKS_pipeline_alert.py**, which loads the archival X-ray data using LoadMasterSources.py and compare it to the PPS file of a new XMM-Newton observation;
3. **StudyMasterSources.py**, which is used for data mining in the archival catalog.
4. **api.py** is the API for the RapidXMM software, used for computation of XMM-Newton upper limits (Ruiz et al. 2021)

This software was developed as part of the XMM2ATHENA project. This project has received funding from the European
Union's Horizon 2020 research and innovation programme under grant agreement n°101004168, the XMM2ATHENA project.

Author: Erwan Quintin, erwan.quintin@irap.omp.eu
