#! /bin/bash

echo "launching stonks catainer"
docker run -d -p 5555:5000 --name stonks --restart unless-stopped  \
     	   --mount source=/data-local/docker_services/stonks/Data,target=/home/stonks/Data,type=bind  \
	   --mount source=/data-local/docker_services/stonks/sessions,target=/home/stonks/sessions,type=bind \
	   -t stonks
echo "launched"
