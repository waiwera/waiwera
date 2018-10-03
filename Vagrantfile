# -*- mode: ruby -*-
# vi: set ft=ruby :

VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.box = "bento/debian-9"
  # If you run into issues with Ansible complaining about executable permissions,
  # comment the following statement and uncomment the next one.
  config.vm.synced_folder ".", "/vagrant/"
  # config.vm.synced_folder ".", "/vagrant", mount_options: ["dmode=700,fmode=600"]
  config.vm.define :waiwera, primary: true do |waiwera|
    waiwera.vm.network :forwarded_port, host: 8088, guest: 8080
    waiwera.vm.network :forwarded_port, host: 5050, guest: 5000
    waiwera.vm.network :forwarded_port, host: 2201, guest: 22, id: "ssh", auto_correct: true
    waiwera.vm.network "private_network", ip: "192.168.50.91"
    waiwera.vm.provision :shell, path: "bootstrap.sh"
    waiwera.vm.provision :shell, path: "waiwera-setup.sh", privileged: false
    waiwera.vm.provision :shell, path: "fruit_processor.sh"
    waiwera.vm.provision :shell, path: "waiwera-testing.sh", privileged: false
    waiwera.vm.hostname = "waiwera"
    waiwera.vm.provision "ansible_local" do |ansible|
      ansible.playbook = "ansible/main.yml"
    end
  end
  config.vm.provider :virtualbox do |v|
    v.gui = true
    v.memory = 4096
  end
  if Vagrant.has_plugin?("vagrant-cachier")
    config.cache.scope = :box
  end
end
