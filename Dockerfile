FROM phusion/baseimage:0.9.15

RUN apt-get update && \
    apt-get install -y python python-dev python-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install ansible
RUN ansible-galaxy install \
    Ansibles.hostname \
    Ansibles.apt \
    Ansibles.build-essential \
    Ansibles.perl \
    Ansibles.monit \
    ANXS.nginx
ADD site.yml /srv/ansible/site.yml

CMD ["/sbin/my_init"]
