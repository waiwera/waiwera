echo "Installing ssh..."
apt install ssh

echo "Installing Ansible..."
apt update
apt upgrade
apt install -y software-properties-common
apt-add-repository ppa:ansible/ansible
apt update
apt install -y ansible
cp /vagrant/ansible/ansible.cfg /etc/ansible/ansible.cfg

echo "Installing Waiwera requirements..."

apt install -y gfortran libopenmpi-dev bison flex cmake git valgrind cvs zlib1g-dev openmpi-bin python g++ vim screen cvs gem rake
# long  build
# apt install -y  xfce4 virtualbox-guest-dkms virtualbox-guest-utils virtualbox-guest-x11

# Python stuff
apt install -y python-pip  python-virtualenv python-h5py virtualenvwrapper python-numpy python-shapely python-future python-pil python-matplotlib python-docutils python-scipy

# pip
# pip install meshio fruitpy

# Enable xfce4 without sudo
echo 'allowed_users=anybody' >> /etc/X11/Xwrapper.config
