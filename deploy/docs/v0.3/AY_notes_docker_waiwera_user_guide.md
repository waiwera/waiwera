# install docker

    https://docs.docker.com/docker-for-windows/install/

On Windows 7 use Docker Toolbox (this one will also use/install VirtualBox VM), Windows 10 install newer docker.  Please read the first few paragraph on the previous link.

# unzip and load waiwera image (only once)

Windows/Mac/Linux:

    docker load -i waiwera-phusion-debian-0.3.image

Linux/Mac?

    gzip -cd waiwera-phusion-debian-0.3.tar.gz | docker load    

# adjust resources on pre-win10 docker (only once)

    docker-machine stop
    /c/Program\ Files/Oracle/VirtualBox/VBoxManage modifyvm default --cpus 12
    /c/Program\ Files/Oracle/VirtualBox/VBoxManage modifyvm default --memory 10240

Also share drive, wherever you possibly wants to run waiwera

    docker-machine start
    docker-machine env

# use one of the following script to run model, in corresponding commandline environment

It is recommended to copy these to somewhere you have PATH set to.
    
    waiwera.sh
    waiwera.bat
    waiwera.ps1

Run by issuing the following command, anywhere you like (as long as the drives are shared)

    waiwera input_file.json 4

# to kill a job cleanly is a slightly trickier

After Ctrl-C, run:

    docker ps

Find the CONTAINER ID of the one that you ran and kill it, for example:

    docker kill 23af998dc4

