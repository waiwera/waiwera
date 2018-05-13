# Geothermal VM Base Setup

This is the base level setup needed to start deploying software to a server

Vagrant needs to be installed for this to work. See step below to setup vagrant:

1. Install VirtualBox
1. Install Vagrant

Run `vagrant up` and a vm called dm will be brought up to run ansible off

- `vagrant up` - starts all machines in the configuration, `vagrant up dm` starts a machine called dm
- `vagrant ssh dm` - connects via ssh to the machine
- `vagrant suspend dm` - suspends the machine called dm

Current machine:

- `dm` - has a basic ansible configuration

ssh_key files are encrypted using ansible vault
