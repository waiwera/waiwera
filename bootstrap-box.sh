echo "Installing ssh..."
apt install ssh

echo "Installing Ansible..."
apt-get update
apt-get upgrade
apt-get install -y software-properties-common
apt-add-repository ppa:ansible/ansible
apt-get update
apt-get install -y ansible
cp /vagrant/ansible/ansible.cfg /etc/ansible/ansible.cfg

echo "Installing Waiwera requirements..."

sudo apt-get install -y gfortran libopenmpi-dev bison flex cmake git valgrind cvs zlib1g-dev openmpi-bin python g++ vim screen cvs gem rake virtualbox-guest-dkms virtualbox-guest-utils virtualbox-guest-x11

# Python stuff
apt-get install -y python-pip  python-virtualenv python-h5py virtualenvwrapper python-numpy python-shapely python-future python-pil python-matplotlib python-docutils python-scipy

pip install meshio
