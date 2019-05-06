$p = "/",(($pwd -replace "\\","/") -replace ":","").ToLower().Trim("/") -join ""
# echo $p $args[0] $args[1]
docker run -v ${p}:/data -w /data waiwera-phusion-debian mpiexec -np $args[1] /home/mpirun/waiwera/dist/waiwera $args[0]
