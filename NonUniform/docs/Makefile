IMAGE := komainu

help:
	@echo "help        : show this message"
	@echo "build-image : build docker image"
	@echo "check       : clippy and rustfmt check"
	@echo "run         : run application to render document"

build-image:
	docker build -t $(IMAGE) .

check:
	docker run \
    --rm \
    -v ${PWD}/Cargo.lock:/app/Cargo.lock:ro \
    -v ${PWD}/Cargo.toml:/app/Cargo.toml:ro \
    -v ${PWD}/src:/app/src \
    $(IMAGE) \
    sh -c "cargo clippy --all-targets && cargo fmt"

run:
	mkdir -p output
	docker run \
    --rm \
    -v ${PWD}/Cargo.lock:/app/Cargo.lock:ro \
    -v ${PWD}/Cargo.toml:/app/Cargo.toml:ro \
    -v ${PWD}/src:/app/src:ro \
    -v ${PWD}/template:/app/template:ro \
    -v ${PWD}/output:/app/output \
    $(IMAGE) \
    cargo run --release

.PHONY: help build-image check run
