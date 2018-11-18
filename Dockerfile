FROM phusion/baseimage:0.11

RUN apt update && \
    apt install -y python python-dev python-pip && \
    apt clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install ansible
ADD ansible /srv/ansible
RUN /usr/local/bin/ansible-playbook /srv/ansible/docker.yml

CMD ["/sbin/my_init"]
