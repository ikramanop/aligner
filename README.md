# Aligner: A tool for sequence alignment

#### Developer: Lev Marder

#### NRNU MEPhI, Spring-Summer 2021

Supplied in form of cli-application and web-service (dispatcher + nodes)

### Building

Use `make build-cli-docker` for building docker image with aligner-cli.

Use `docker run -v $(pwd):/data aligner-cli:$branch $args` to run tool after building.

Use `--help` option after build to see how to work with the tool.
