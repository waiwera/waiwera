# Geothermal VM Base Setup

This is the base level setup needed to start deploying software to a server

Vagrant needs to be installed for this to work. See step below to setup vagrant:

1. Install VirtualBox from: https://www.virtualbox.org/wiki/Downloads
1. Install Vagrant from: https://www.vagrantup.com/downloads.html
1. Download the vagrant setup from: link to be added

Run `vagrant init` and within the geo-deploy directory

- `vagrant up` - starts all machines in the configuration, `vagrant up waiwera` starts a machine called `waiwera`
- `vagrant ssh waiwera` - connects via ssh to the machine
- `cd ~/supermodels-test/`
- `python unit_tests.py`
- `vagrant suspend waiwera` - suspends the machine called dm

Current machine:

- `waiwera` - has a basic ansible configuration

ssh_key files are encrypted using ansible vault
