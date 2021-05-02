export CLI_NAME=aligner-cli
export BRANCH=$(shell git symbolic-ref --short HEAD)

build-cli-docker:
	docker build -t ${CLI_NAME}:${BRANCH} -f Dockerfile-cli .