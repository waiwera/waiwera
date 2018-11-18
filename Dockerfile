FROM phusion/baseimage:0.11
# FROM debian/9

RUN apt update && \
    apt install -y python python-dev python-pip && \
    apt clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install ansible
ADD ansible /srv/ansible
RUN /usr/local/bin/ansible-playbook /srv/ansible/main.yml -v
