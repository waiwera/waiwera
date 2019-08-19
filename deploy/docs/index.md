# Waiwera Deployment

For full documentation visit [mkdocs.org](https://mkdocs.org).

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs help` - Print this help message.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.

# Waiwera box setup

Setup a virtualbox vm for exporting to vagrant

* Install Ubuntu 18.04 Desktop
  * Username: `vagrant`
  * Password: `vagrant`
  * network: `bridged`


Bootstrap the box
```
sudo apt install ssh
```
```
ssh vagrant@waiwera.local
scp bootstrap.sh vagrant@waiwera.local:
```

```
cd ~
sudo sh bootstrap.sh
```

Shutdown VM and package it up using Vagrant

start xfce4 using these instructions
https://stackoverflow.com/questions/18878117/using-vagrant-to-run-virtual-machines-with-desktop-environment
