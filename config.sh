#!/bin/sh
cd build
CC=mpicc FC=mpif90 cmake ..
cd ..
