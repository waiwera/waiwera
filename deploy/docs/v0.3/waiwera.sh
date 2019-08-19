#!/bin/bash
echo running $1 in $(pwd)
# docker run -v $(pwd):/data -w /data alpine ls
docker run -v $(pwd):/data -w /data waiwera-phusion-debian mpiexec -np $2 /home/mpirun/waiwera/dist/waiwera $1
