### CREDO2, required for benchmark testing
cd ~
git clone https://angusyeh@bitbucket.org/uoa-geo/credo2.git
# need access
cd credo2
git checkout refactor-multi-simulators


### PyTOUGH required by benchmark testing
cd ~
git clone https://github.com/acroucher/PyTOUGH.git


### AUTOUGH2, required for benchmark testing
cd ~
git clone https://angusyeh@bitbucket.org/uoa-geo/autough2.git
# need access
cd autough2
git checkout 2.42
make -f makefile_gcc_linux
cp AUTOUGH2_42D ~/bin/
