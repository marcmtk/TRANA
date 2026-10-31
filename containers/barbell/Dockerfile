# barbell built from the marcmtk fork, which adds the SQK-MAB114-24 16S/ITS presets
FROM rust:1.90-slim-bookworm AS build
ARG BARBELL_COMMIT=086ab1dea26d5d3ba572fccf2fe7b6f435ec6657
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates build-essential
WORKDIR /src
RUN git init -q . \
    && git remote add origin https://github.com/marcmtk/barbell.git \
    && git fetch -q --depth 1 origin "${BARBELL_COMMIT}" \
    && git checkout -q FETCH_HEAD
# The repo defaults to target-cpu=native, which would not be portable between hosts
RUN RUSTFLAGS="-C target-cpu=x86-64-v3" cargo build --release --locked --bin barbell \
    && cp target/release/barbell /usr/local/bin/barbell

FROM debian:bookworm-slim
ARG BARBELL_COMMIT=086ab1dea26d5d3ba572fccf2fe7b6f435ec6657
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
COPY --from=build /usr/local/bin/barbell /usr/local/bin/barbell
ENV BARBELL_COMMIT=${BARBELL_COMMIT}
LABEL org.opencontainers.image.source="https://github.com/marcmtk/barbell" \
      org.opencontainers.image.revision="${BARBELL_COMMIT}"
