# -*- mode: ruby -*-
# vi: set ft=ruby :

VAGRANTFILE_API_VERSION = "2"
Vagrant.require_version ">= 2.1.4"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.box = "bento/debian-9"
  config.vm.synced_folder ".", "/vagrant/"
  # Problems with meson locale and UTF-8
  ENV['LC_ALL']="en_US.UTF-8"
  config.vm.define :waiwera, primary: true do |waiwera|
    waiwera.vm.network :forwarded_port, host: 8088, guest: 8080
    waiwera.vm.network :forwarded_port, host: 5050, guest: 5000
    waiwera.vm.network :forwarded_port, host: 2201, guest: 22, id: "ssh", auto_correct: true
    waiwera.vm.network "private_network", ip: "192.168.50.91"
    # waiwera.vm.provision :shell, path: "bootstrap.sh"
    # waiwera.vm.provision :shell, path: "waiwera-setup.sh", privileged: false
    # waiwera.vm.provision :shell, path: "fruit_processor.sh"
    # waiwera.vm.provision :shell, path: "waiwera-testing.sh", privileged: false
    waiwera.vm.provision "ansible" do |ansible|
      ansible.limit = 'waiwera'
      ansible.verbose = "vv"
      ansible.playbook = "ansible/packer.yml"
      ansible.raw_arguments = [
        '-e "waiwera_user=' + ENV['W_USER'] + '"',
        '-e "waiwera_pwd='+ ENV['W_PWD'] +'"',
        '--skip-tags=docker,clean',
      ]
    end
    # config.vm.provision "docker" do |d|
    #   d.build_image "/vagrant/"
    # end
    waiwera.vm.hostname = "waiwera"
  end
  config.vm.provider :virtualbox do |v|
#    v.gui = true
    v.memory = 4096
  end
  if Vagrant.has_plugin?("vagrant-cachier")
    config.cache.scope = :box
  end
end
