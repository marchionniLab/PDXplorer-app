default:
  @just --choose;

run-local:
  #! /usr/bin/env -S bash -i
  \R -q -s -e \
    'shiny::runApp(host="127.0.0.1",port=8081,launch.browser=FALSE)';

sync server='omega':
  #! /usr/bin/env -S bash -i
  rsync -az -e 'ssh' -r --delete --info=progress2 \
  --exclude 'renv' --exclude '.git' \
  {{justfile_directory()}} {{server}}:projects;

build-container:
  #! /usr/bin/env -S bash -i
  docker buildx build -t pdxplorer-app:latest .;

build-container-mac:
  #! /usr/bin/env -S bash -i
  docker buildx build --platform=linux/amd64 -t pdxplorer-app:latest .;

run-container port='8081':
  #! /usr/bin/env -S bash -i
  docker stop pdxplorer-app;
  docker rm pdxplorer-app;
  docker run -d \
    --restart always \
    -p {{port}}:80 \
    --name pdxplorer-app \
    pdxplorer-app:latest;

monitor:
  #! /usr/bin/env -S bash -i
  docker logs -f pdxplorer-app;
