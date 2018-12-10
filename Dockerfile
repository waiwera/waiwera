FROM debian:9

ARG waiwera_user
ARG waiwera_pwd

ENV PETSC_DIR=/opt/app/petsc
ENV PETSC_ARCH=arch-linux2-c-debug
ENV LD_LIBRARY_PATH="/opt/app/lib"
ENV PYTHONPATH=/"/opt/app/PyTOUGH:/opt/app/credo2:/opt/app/supermodels-test/utils"
ENV PATH="$PATH:/opt/app/bin:/opt/app/supermodels-test/dist"

RUN apt update && \
    apt install -y python python-dev python-pip && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*  /root/.cache/*

RUN pip install pip
RUN pip install ansible

ADD ansible /srv/ansible
WORKDIR /srv/

RUN /usr/local/bin/ansible-playbook -c local ansible/site.yml -e "waiwera_user=${waiwera_user} waiwera_pwd=${waiwera_pwd}"

WORKDIR /opt/app/supermodels-test

RUN pip uninstall ansible -y && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/* /srv/ansible /usr/share/doc/*
RUN for dep in $(pip show ansible | grep Requires | sed 's/Requires: //g; s/,//g'); do pip uninstall -y $dep; done

RUN find / -uid|-gid 165586 -ls
