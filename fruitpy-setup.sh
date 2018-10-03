### FRUITpy, required for unit testing
cd ~
git clone https://github.com/acroucher/FRUITPy.git fruitpy
sudo pip install fruitpy/


### FRUIT, required for unit testing
cd ~
cvs -d:pserver:anonymous@fortranxunit.cvs.sourceforge.net:/cvsroot/fortranxunit login
# password just empty
cvs -z3 -d:pserver:anonymous@fortranxunit.cvs.sourceforge.net:/cvsroot/fortranxunit co -P fruit
cd fruit/fruit_processor_gem
rake install


### FRUITpy, required for unit testing
cd ~
cp fruitpy/fruit_makefile fruit/makefile
cd fruit
make
make install
cp *.mod ~/include
