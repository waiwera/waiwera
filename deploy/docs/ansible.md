# Ansible
Things to use with Ansible

```
# If you run into issues with Ansible complaining about executable permissions,
# comment the following statement and uncomment the next one.

# config.vm.synced_folder ".", "/vagrant", mount_options: ["dmode=700,fmode=600"]
```

Within   `config.vm.define :waiwera, primary: true do |waiwera|`
```
# waiwera.vm.provision "ansible_local" do |ansible|
#   ansible.playbook = "ansible/main.yml"
# end
```
