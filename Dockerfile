FROM python:3.8-slim
LABEL David Corre <corre@iap.fr>

RUN  \
   apt-get update \
   && apt-get install -y build-essential openssh-client git curl default-jdk
RUN apt-get install -y libtk-img-dev 
RUN rm -rf /var/lib/apt/lists/* 

# User setup
ARG svom_gid=998
ARG svom_uid=208
ENV USR svom
ENV SVOM_GID $svom_gid
ENV SVOM_UID $svom_uid
RUN addgroup -gid $SVOM_GID $USR
RUN adduser --system --uid $SVOM_UID --gid $SVOM_GID $USR
# Now logon as svom user
USER $USR:$USR
WORKDIR /home/$USR

ADD ./requirement.txt .

ENV PYTHONPATH /home/svom/python:$PYTHONPATH
ENV PATH /home/svom/.local/bin/:$PATH

RUN  \
   python -m pip install --upgrade pip \
   && python -m pip install -U -r requirement.txt \
   && rm -fr ${HOME}/.cache/pip 

ADD ./python python
ADD ./stilts stilts
ADD ./doc doc
RUN mkdir Data
# Run the production server
CMD gunicorn --bind 0.0.0.0:5000 rest_api.wsgi_api:application
#CMD bash
