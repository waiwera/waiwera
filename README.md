# Geothermal VM Base Setup

This is the base level setup needed to start deploying software to a server

## Local Machine Setup
Git, Vagrant, and VirtualBox need to be installed for this to work. The repository which sets up the test environment also needs to be cloned. See step below to setup vagrant:

1. Install Git from: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
1. Install VirtualBox from: https://www.virtualbox.org/wiki/Downloads
1. Install Vagrant from: https://www.vagrantup.com/downloads.html (Make sure version >= 2)
1. Download the vagrant setup from:
```
git clone https://timharton@bitbucket.org/timharton/geo-deploy.git
```

**Make sure to use http not https**

## Vagrant Virtual Machine setup

Open the command line for your operating system and follow these steps:

- Navigate to the he directory geo-deploy was cloned into and run `vagrant init`
- This only has to be done once to tell the system that this is vagrant deployment directory
- `vagrant up` - starts a VM (virtual machine) named `waiwera`
- This could take some time as it builds and installs all the dependencies required
- `vagrant ssh` - connects via ssh to the machine. Within the `ssh` session run:
- `cd ~/supermodels-test/`
- `python unit_tests.py`, all tests should pass
- `exit` will leave the ssh connection
- You are now done you can either destroy or suspend the `waiwera` VM
- `vagrant suspend` - suspends the VM
- `vagrant destroy` - destroys the VM

All commands should be run be run from the deployment directory. Otherwise the name of the machine can be added to the end of the vagrant command i.e. `vagrant suspend waiwera`.
