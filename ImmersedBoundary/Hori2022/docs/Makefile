IMAGE := naokihori/komainu:latest

# Source is mounted read-only and the container has no network access (no dependencies).
DOCKER_RUN := docker run --rm --network none

help:
	@echo "help        : show this message"
	@echo "build-image : build docker image"
	@echo "check       : clippy and rustfmt check"
	@echo "run         : run application to render document"

build-image:
	docker build -t $(IMAGE) .

check:
	$(DOCKER_RUN) -v ${PWD}:/app $(IMAGE) sh -c "cargo clippy --all-targets && cargo fmt"

run:
	$(DOCKER_RUN) -v ${PWD}:/app $(IMAGE) cargo run --release

.PHONY: help build-image check run
