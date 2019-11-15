
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

```
docker run --rm --group-add audio --group-add nogroup --group-add 777 busybox id
uid=0(root) gid=0(root) groups=10(wheel),29(audio),99(nogroup),777

docker container run --rm -it \
  -v $(app):/app \                          # Mount the source code
  --workdir /app \                          # Set the working dir
  --user 1000:1000 \                        # Run as the given user
  my-docker/my-build-environment:latest \   # Our build env image
  make assets
 ```

````
docker container run --rm -it \
  --user 1000:1000 \                        # Run as the given user
  thar102/waiwera:latest \   # Our build env image
  python /opt/app/supermodels-test unit_tests.py
```
