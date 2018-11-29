# FROM phusion/baseimage:0.11
FROM debian:9
RUN groupadd -r waiwera && useradd --no-log-init -r -g waiwera waiwera

RUN apt update && \
    apt install -y python python-dev python-pip && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*  /root/.cache/*

RUN pip install pip
RUN pip install ansible

ADD ansible /srv/ansible
RUN /usr/local/bin/ansible-playbook -c local /srv/ansible/site.yml  -v

RUN pip uninstall ansible -y && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/* /srv/ansible /usr/share/doc/*
RUN for dep in $(pip show ansible | grep Requires | sed 's/Requires: //g; s/,//g'); do pip uninstall -y $dep; done
