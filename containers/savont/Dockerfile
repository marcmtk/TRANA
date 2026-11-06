# savont built from the marcmtk fork, which adds the SBDI curated GTDB database
FROM rust:1.90-slim-bookworm AS build
ARG SAVONT_COMMIT=f828ffd4bcc84f7d4fa1744be9803fd392512d48
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates build-essential clang libclang-dev cmake pkg-config zlib1g-dev
WORKDIR /src
RUN git init -q . \
    && git remote add origin https://github.com/marcmtk/savont.git \
    && git fetch -q --depth 1 origin "${SAVONT_COMMIT}" \
    && git checkout -q FETCH_HEAD
RUN RUSTFLAGS="-C target-cpu=x86-64-v3" cargo build --release --locked --bin savont \
    && cp target/release/savont /usr/local/bin/savont

FROM debian:bookworm-slim
ARG SAVONT_COMMIT=f828ffd4bcc84f7d4fa1744be9803fd392512d48
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
COPY --from=build /usr/local/bin/savont /usr/local/bin/savont
ENV SAVONT_COMMIT=${SAVONT_COMMIT}
LABEL org.opencontainers.image.source="https://github.com/marcmtk/savont" \
      org.opencontainers.image.revision="${SAVONT_COMMIT}"
