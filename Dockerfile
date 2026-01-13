# Multi-stage build for smaller image
FROM python:3.11-slim as builder

WORKDIR /build

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy and install package
COPY . .
RUN pip wheel --no-cache-dir --wheel-dir /wheels -e .


FROM python:3.11-slim

LABEL maintainer="BSAseq Authors <bsaseq@example.com>"
LABEL description="bsaseq: Bulk Segregant Analysis for QTL mapping"
LABEL version="1.0.0"
LABEL org.opencontainers.image.source="https://github.com/username/bsaseq"
LABEL org.opencontainers.image.documentation="https://bsaseq.readthedocs.io"
LABEL org.opencontainers.image.licenses="MIT"

WORKDIR /data

# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# Copy wheels and install
COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir /wheels/* && rm -rf /wheels

# Install snpEff (optional, makes image larger)
ARG INSTALL_SNPEFF=false
RUN if [ "$INSTALL_SNPEFF" = "true" ]; then \
    apt-get update && apt-get install -y --no-install-recommends \
    default-jre-headless \
    wget \
    unzip \
    && wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
    && unzip snpEff_latest_core.zip -d /opt \
    && rm snpEff_latest_core.zip \
    && ln -s /opt/snpEff/exec/snpEff /usr/local/bin/snpEff \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*; \
    fi

# Create non-root user
RUN useradd -m -s /bin/bash bsaseq
USER bsaseq

ENTRYPOINT ["bsaseq"]
CMD ["--help"]
