#
# STONK server image
#
# ABSOLUTEPATHTODATA/Data subfolders must all have a write permission
# docker run -p 5600:5000 --mount source=ABSOLUTEPATHTODATA/Data,target=/home/stonks/Data,type=bind -it stonks
#
#
FROM python:3.8-slim
LABEL Laurent Michel <laurent.michel@astro.unistra.fr>

# install Linux packages
# java takes a piece of time
RUN  \
   apt-get update \
   && apt-get install -y build-essential openssh-client git curl default-jdk
RUN apt-get install -y libtk-img-dev 
RUN rm -rf /var/lib/apt/lists/* 

# User setup
ARG stonks_gid=998
ARG stonks_uid=208
ENV USR stonks
ENV stonks_GID $stonks_gid
ENV stonks_UID $stonks_uid
RUN addgroup -gid $stonks_GID $USR
RUN adduser --system --uid $stonks_UID --gid $stonks_GID $USR
# Now logon as stonks user
USER $USR:$USR
WORKDIR /home/$USR

# get the list of required Python packages
ADD ./requirement.txt .

ENV PYTHONPATH /home/stonks/python:$PYTHONPATH
ENV PATH /home/stonks/.local/bin/:$PATH

# install Python packages
RUN  \
   python -m pip install --upgrade pip \
   && python -m pip install -U -r requirement.txt \
   && rm -fr ${HOME}/.cache/pip 

# Get Python sources and STILTS Java tool 
ADD ./python python
ADD ./stilts stilts
# Get tyhe doc static files
ADD ./doc doc
ADD ./static static

# Set the mount point for the science stuff
RUN mkdir Data
RUN mkdir sessions

# Run the production server
CMD gunicorn --bind 0.0.0.0:5000 --timeout 600 rest_api.wsgi_api:application

