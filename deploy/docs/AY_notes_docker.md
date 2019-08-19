# Working on 16/10/2018 1:32:34 p.m.

I had difficulties getting MPI to work when I directly use ubuntu or debian, so I started looking into those "base images" available online for MPI applications.  I had tried:

https://github.com/NLKNguyen/alpine-mpich
https://github.com/oweidner/docker.openmpi
https://github.com/phusion/baseimage-docker

The oweidner/docker.openmpi one runs the example in a strange way (one CPU process in each container, needs a head and several nodes to run MPI), and I couldn't get petsc to work even with `make check`.  The NLKNguyen/alpine-mpich one was not able to configure Petsc properly, which  failed at making pnetcdf.

I was only able to get the phusion one to work.

* phusion-base\
* petsc-phusion-debian\
* waiwera-phusion-debian\

The phusion baseimage is based on ubuntu, and its gfortran compiler (5.4.0) does NOT work well with waiwera's list module (unit tests fail with list tags, hence lots dependent modules also fails).  So I had to recompile phusion's base image from Debian (stretch has gfortran 6.3.0)

Rebuild the base image following the instructions: https://github.com/phusion/baseimage-docker#building

I had to modify a few things in there to get it to run:
```
    en-354404+cyeh015@en-354404 /cygdrive/d/_software/docker/phusion-base
    $ git diff
    diff --git a/image/Dockerfile b/image/Dockerfile
    index da5b41d..45ed558 100644
    --- a/image/Dockerfile
    +++ b/image/Dockerfile
    @@ -1,4 +1,4 @@
    -FROM ubuntu:18.04
    +FROM debian:stretch
     MAINTAINER Phusion <info@phusion.nl>

     COPY . /bd_build
    diff --git a/image/prepare.sh b/image/prepare.sh
    index c2926a5..99ff479 100755
    --- a/image/prepare.sh
    +++ b/image/prepare.sh
    @@ -41,8 +41,9 @@ $minimal_apt_get_install software-properties-common
     apt-get dist-upgrade -y --no-install-recommends -o Dpkg::Options::="--force-confold"

     ## Fix locale.
    -$minimal_apt_get_install language-pack-en
    +$minimal_apt_get_install locales locales-all
     locale-gen en_US
    +locale-gen en_US.UTF-8
     update-locale LANG=en_US.UTF-8 LC_CTYPE=en_US.UTF-8
     echo -n en_US.UTF-8 > /etc/container_environment/LANG
     echo -n en_US.UTF-8 > /etc/container_environment/LC_CTYPE
```
Then after building image, I had to copy it out of the virtualbox:
```
    vagrant@ubuntu-14:/vagrant$ docker save phusion-baseimage-debian > phusion-baseimage-debian
```
Then load it back on my normal docker environment registry (windows):
```
    cat phusion-baseimage-debian | docker load
```
Now I have a new baseimage based on debian, I then build petsc and waiwera using the customised image, see `petsc-phusion-debian\` and `waiwera-phusion-debian\`.

Now all unit tests pass.  Run the image with command:
    docker run -w="/home/mpirun/test_medium" waiwera-phusion-debian mpiexec -np 8 ../waiwera/dist/waiwera 2DM002.json

And it seems to work fine, with multiple CPU using up proper resources.
