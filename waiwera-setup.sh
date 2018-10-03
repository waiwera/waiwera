echo export PETSC_DIR=/home/vagrant/petsc >> /home/vagrant/.bashrc
echo export PETSC_ARCH=arch-linux2-c-debug >> /home/vagrant/.bashrc
echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vagrant/lib >> /home/vagrant/.bashrc
echo export PYTHONPATH=/home/vagrant/PyTOUGH:/home/vagrant/credo2:/home/vagrant/supermodels-test/utils >> /home/vagrant/.bashrc
echo export PATH=$PATH:/home/vagrant/bin:/home/vagrant/supermodels-test/dist >> /home/vagrant/.bashrc
source /home/vagrant/.bashrc

pip install meshio

### PETSc, required by Waiwera
git clone https://angusyeh@bitbucket.org/petsc/petsc.git /home/vagrant/petsc
cd /home/vagrant/petsc
# git checkout next
git checkout master
# ./configure --download-fblaslapack --download-triangle --download-netcdf --download-exodusii --download-pnetcdf --download-hdf5 --download-ptscotch --download-chaco
./configure --download-fblaslapack --download-triangle --download-netcdf --download-exodusii --download-pnetcdf --download-hdf5 --download-ptscotch --download-chaco --with-zlib
make

# --useThreads=0
# --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-shared-libraries

mkdir /home/vagrant/bin
mkdir /home/vagrant/lib
mkdir /home/vagrant/include

### FRUITpy, required for unit testing
git clone https://github.com/acroucher/FRUITPy.git /home/vagrant/fruitpy
pip install /home/vagrant/fruitpy/

### FRUIT, required for unit testing
cd /home/vagrant
cvs -d:pserver:anonymous@fortranxunit.cvs.sourceforge.net:/cvsroot/fortranxunit login
# password just empty
cvs -z3 -d:pserver:anonymous@fortranxunit.cvs.sourceforge.net:/cvsroot/fortranxunit co -P fruit
