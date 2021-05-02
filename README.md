# Aligner: A tool for sequence alignment

#### Developer: Lev Marder

#### NRNU MEPhI, Spring 2021

Currently supplied in cli only.

### Building

Use `make build-cli-docker` for building docker image with aligner-cli.

Use `docker run -v $(pwd):/data aligner-cli:$branch $args` to run tool after building.

Use `--help` option after build to see how to work with the tool.
