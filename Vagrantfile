# -*- mode: ruby -*-
# vi: set ft=ruby :

VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.box = "ubuntu/xenial64"
  # If you run into issues with Ansible complaining about executable permissions,
  # comment the following statement and uncomment the next one.
  config.vm.synced_folder ".", "/vagrant"
  # config.vm.synced_folder ".", "/vagrant", mount_options: ["dmode=700,fmode=600"]
  config.vm.provider "virtualbox" do |v|
    v.memory = 512
  end
  config.vm.define :dm, primary: true do |dm|
    dm.vm.network :forwarded_port, host: 8088, guest: 8080
    dm.vm.network :forwarded_port, host: 5050, guest: 5000
    dm.vm.network :forwarded_port, host: 2201, guest: 22, id: "ssh", auto_correct: true
    dm.vm.network "private_network", ip: "192.168.50.91"
    dm.vm.provision "shell", path: "bootstrap.sh"
    dm.vm.provision :shell, inline: 'ansible-playbook /vagrant/ansible/dm.yml -c local -v'
    dm.vm.hostname = "dm"
  end
  if Vagrant.has_plugin?("vagrant-cachier")
    config.cache.scope = :box
  end
end