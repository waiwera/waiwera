
List images and run image after build

```
docker images
docker run {{ image_id }}
```

List running containers and get a command line interface

```
docker ps
docker exec -it d6529ce008f3 /bin/bash
```

```
docker run -i -t dc0a1991fa15 /bin/bash
```

build waiwera container
```
docker build -t waiwera /vagrant/
```
