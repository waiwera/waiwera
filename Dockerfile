#FROM phusion/baseimage
FROM debian:bookworm-slim

USER root

# ARG git_user
# ARG git_pwd
ARG base_dir=/opt

RUN DEBIAN_FRONTEND=noninteractive \
    apt-get update && \
    apt-get install -y \
        locales \
        wget && \
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

# install gosu for a better su+exec command
ARG GOSU_VERSION=1.17
RUN dpkgArch="$(dpkg --print-architecture | awk -F- '{ print $NF }')" \
 && wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$dpkgArch" \
 && chmod +x /usr/local/bin/gosu \
 && gosu nobody true

#ADD install/ansible /ansible
ADD . /opt/waiwera
#RUN rm -rf /opt/waiwera/external && rm -rf /opt/waiwera/build

RUN DEBIAN_FRONTEND=noninteractive \
    apt-get update && \
    apt-get install -y --no-install-recommends  \
        apt-utils \
        python3-minimal \
        python3-venv && \
    python3 -m venv ~/env && \
    . ~/env/bin/activate && \
    ~/env/bin/pip install \
        setuptools \
        wheel \
        ansible && \
    ~/env/bin/ansible-playbook \
        --connection=local \
        /opt/waiwera/install/ansible/docker.yml \
        -v && \
    deactivate && \
    rm -rf ~/env && \
    apt-get remove -y python3-venv && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    apt-get clean -y && \
    rm -rf \
        /var/lib/apt/lists/* \
        /tmp/* \
        /var/tmp/* \
        /usr/share/doc/* \
        /root/.pip/cache/*

#RUN rm -r /ansible

# USER 1000:1000
# Setup docker group
#RUN groupadd -r docker \
# && usermod -aG docker waiwera

#COPY entrypoint.sh ${base_dir}/waiwera/

WORKDIR ${base_dir}/waiwera
RUN mkdir /data && chown waiwera:waiwera /data
#ARG entry=${base_dir}/waiwera/entrypoint.sh
#ENV entry $entry

# ARGs can't be used with ENTRYPOINT command
ENTRYPOINT [ "/opt/waiwera/entrypoint.sh" ]
