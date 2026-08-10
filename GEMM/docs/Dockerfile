FROM rust:latest

RUN rustup component add clippy rustfmt

WORKDIR /app

ENV CARGO_NET_GIT_FETCH_WITH_CLI=true

RUN echo '#!/bin/sh\n\
if [ -n "$GH_TOKEN" ]; then\n\
  git config --global url."https://x-access-token:${GH_TOKEN}@github.com/".insteadOf "https://github.com/"\n\
fi\n\
exec "$@"' > /entrypoint.sh && chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
