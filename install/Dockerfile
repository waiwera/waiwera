FROM debian:stretch-slim

# ARG git_user
# ARG git_pwd
ARG base_dir=/opt

RUN DEBIAN_FRONTEND=noninteractive \
    apt-get update && \
    apt-get install -y \
        locales && \
    rm -rf /var/lib/apt/lists/* && \
    localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8 && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf \
        /var/lib/apt/lists/* \
        /tmp/* \
        /var/tmp/* \
        /usr/share/doc/*  \
        /root/.pip/cache/*
ENV LANG en_US.utf8


ADD ansible /ansible

RUN DEBIAN_FRONTEND=noninteractive \
    apt update && \
    apt install -y --no-install-recommends  \
        python-minimal \
        python-pip \
        python-setuptools \
        python-wheel \
        apt-utils && \
    pip install ansible &&\
    usr/local/bin/ansible-playbook \
        --connection=local \
        /ansible/docker.yml \
        -v \
#        -e "base_dir=${base_dir}" \
#        --tags=all,docker \
#        --skip-tags=local,zofu,fson
         && \
    pip uninstall ansible -y && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    apt-get clean -y && \
    rm -rf \
        /var/lib/apt/lists/* \
        /tmp/* \
        /var/tmp/* \
        /usr/share/doc/*  \
        /root/.pip/cache/*

RUN rm -r /ansible

USER 1000:1000

#COPY entrypoint.sh ${base_dir}/waiwera/

WORKDIR ${base_dir}/waiwera
#ARG entry=${base_dir}/waiwera/entrypoint.sh
#ENV entry $entry

# ARGs can't be used with ENTRYPOINT command
ENTRYPOINT [ "/opt/waiwera/entrypoint.sh" ]