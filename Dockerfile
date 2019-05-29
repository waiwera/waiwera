FROM debian:stretch-slim

ARG git_user
ARG git_pwd
ARG app_dir

ENV PETSC_DIR=${app_dir}/waiwera/external/PETSc
ENV PETSC_ARCH=arch-linux2-c-debug
ENV LD_LIBRARY_PATH="${app_dir}/lib"
ENV PYTHONPATH="${app_dir}/PyTOUGH:${app_dir}/credo2:${app_dir}/waiwera/utils"
ENV PATH="$PATH:${app_dir}/bin:${app_dir}/waiwera/dist"
ENV PKG_CONFIG_PATH="$PKG_CONFIG_PATH:${app_dir}/lib/pkgconfig"

RUN apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* && \
    localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8 && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*  /root/.pip/cache/*
ENV LANG en_US.utf8


ADD ansible /ansible

RUN apt update && \
    apt install -y --no-install-recommends python3-minimal python3-pip python3-setuptools python3-wheel && \
    pip3 install ansible &&\
    usr/local/bin/ansible-playbook --connection=local /ansible/local.yml -e "waiwera_user=${git_user}" -e "waiwera_pwd=${git_pwd}" -e  "app_dir=${app_dir}" --skip-tags=vagrant,local,packer,clean -v && \
    pip3 uninstall ansible -y && \
    for dep in $(pip3 show ansible | grep Requires | sed 's/Requires: //g; s/,//g'); do pip3 uninstall -y $dep; done && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*  /root/.pip/cache/*

RUN rm -r /ansible

WORKDIR ${app_dir}/waiwera
