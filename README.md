The Search for Transient Objects in New detections using Known Sources (STONKS) software is meant to allow automatic
comparison between the detections of a new XMM-Newton observation, and all the available archival X-ray data. It will send
out an alert and save the lightcurve and available data for all sources with a long-term variability over a factor of 5.
The updated framework is meant to be used with a single input position and flux state.


This repository contains 3 Python scripts:
1. **LoadSpecificMasterSources.py**, which is used to load the archival X-ray data in the form of a MasterSource object matching the input position;
2. **STONKS_PreComputed_Position_alert.py**, which will compare the new detection to any archival MasterSource, and compute upper-limits if there are no known MasterSource;
3. **api.py** is the API for the RapidXMM software, used for computation of XMM-Newton upper limits (Ruiz et al. 2021)

Further documentation is available in STONKS_Documentation.pdf if needed.

This software was developed as part of the XMM2ATHENA project. This project has received funding from the European
Union's Horizon 2020 research and innovation programme under grant agreement n°101004168, the XMM2ATHENA project.

Author: Erwan Quintin, erwan.quintin@irap.omp.eu

NOTE: In its current form, the software will not work because it will be missing the catalog data. Once a solution has been found to host this massive data in a sustainable fashion, this repository will be updated.

# Launch command

docker run -d -p 5555:5000 --name stonks --restart unless-stopped --mount source=/home/laurent.michel/STONKS/Data,target=/home/stonks/Data,type=bind -t stonks

# Start at boot:
```
[Unit]
Description=STONKS container
After=docker.service
Wants=network-online.target docker.socket
Requires=docker.socket

[Service]
Restart=always
ExecStart=/usr/bin/docker start -a stonks
ExecStop=/usr/bin/docker stop -t 3 stonks

[Install]
WantedBy=multi-user.target
```
