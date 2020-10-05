# Build a development container

Use this Dockerfile to create a container with nightly OpenMS versions.
Depends on the DockerHub builds of [openms/executables](https://hub.docker.com/r/openms/executables).
Build with:

```bash
docker build -t proteomicslfq-dev .
```

Then use the dev profile to use this container.
