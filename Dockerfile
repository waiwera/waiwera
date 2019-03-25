FROM debian:9

ARG waiwera_user
ARG waiwera_pwd
ARG app_dir

ENV PETSC_DIR=${app_dir}/waiwera/external/PETSc
ENV PETSC_ARCH=arch-linux2-c-debug
ENV LD_LIBRARY_PATH="${app_dir}/lib"
ENV PYTHONPATH=/"${app_dir}/PyTOUGH:${app_dir}/credo2:${app_dir}waiwera/utils"
ENV PATH="$PATH:${app_dir}/bin:${app_dir}/waiwera/dist"
ENV PKG_CONFIG_PATH=$"PKG_CONFIG_PATH:${app_dir}/lib/pkgconfig"

RUN apt-get clean && apt-get update && apt-get install -y locales locales-all

ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
RUN locale-gen en_US.UTF-8
RUN dpkg-reconfigure locales

RUN apt update && \
    apt install -y python3-minimal python3-pip && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*  /root/.cache/*

RUN pip3 install pip
RUN pip3 install ansible

ADD ansible /srv/ansible
WORKDIR /srv/

RUN /usr/local/bin/ansible-playbook --connection=local ansible/docker.yml -e "waiwera_user=${waiwera_user}" -e "waiwera_pwd=${waiwera_pwd}" -e "app_dir=${app_dir}" --skip-tags=vagrant

WORKDIR ${app_dir}/waiwera

RUN pip3 uninstall ansible -y && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/* /srv/ansible /usr/share/doc/*

RUN for dep in $(pip3 show ansible | grep Requires | sed 's/Requires: //g; s/,//g'); do pip3 uninstall -y $dep; done
