DOCKER_IMAGE := naokihori/cfd-toolbox-docs:latest
OPTS := --rm -it -v ${PWD}:/app -w /app -e RUST_BACKTRACE=1 --env-file .env --init

build:
	docker run $(OPTS) $(DOCKER_IMAGE) cargo build

clean:
	docker run $(OPTS) $(DOCKER_IMAGE) cargo clean

clippy:
	docker run $(OPTS) $(DOCKER_IMAGE) cargo clippy --tests

fmt:
	docker run $(OPTS) $(DOCKER_IMAGE) cargo fmt

run:
	docker run $(OPTS) $(DOCKER_IMAGE) cargo run

test:
	docker run $(OPTS) $(DOCKER_IMAGE) cargo test -- --nocapture
