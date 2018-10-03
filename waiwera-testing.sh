# export PETSC_DIR=/home/vagrant/petsc
# export PETSC_ARCH=arch-linux2-c-debug
source /root/.bashrc

### FRUITpy, required for unit testing
cp /home/vagrant/fruitpy/fruit_makefile /home/vagrant/fruit/makefile
cd /home/vagrant/fruit
make
make install
cp /home/vagrant/fruit/*.mod /home/vagrant/include


### fson, required by Waiwera
git clone https://github.com/josephalevin/fson.git /home/vagrant/fson
cd /home/vagrant/fson
make
make install

### main simulator waiwera
cd /home/vagrant
git clone https://jposunz:wapiti21@github.com/johnburnell/supermodels-test.git /home/vagrant/supermodels-test
# need access
cd /home/vagrant/supermodels-test
git checkout testing
make
# python unit_tests.py
