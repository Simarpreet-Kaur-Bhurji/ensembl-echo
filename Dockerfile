FROM python:3.10-slim

WORKDIR /app

# System dependencies needed by ete3 and matplotlib
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    gcc \
    libgl1 \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install MMseqs2 — pick the right binary for the target CPU architecture
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "arm64" ]; then \
      wget -q https://mmseqs.com/latest/mmseqs-linux-arm64.tar.gz \
      && tar xzf mmseqs-linux-arm64.tar.gz \
      && mv mmseqs/bin/mmseqs /usr/local/bin/mmseqs \
      && rm -rf mmseqs mmseqs-linux-arm64.tar.gz; \
    else \
      wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz \
      && tar xzf mmseqs-linux-avx2.tar.gz \
      && mv mmseqs/bin/mmseqs /usr/local/bin/mmseqs \
      && rm -rf mmseqs mmseqs-linux-avx2.tar.gz; \
    fi

# Verify installations
RUN python -c "import pandas, numpy, matplotlib, duckdb, ete3; print('All packages OK')"
RUN mmseqs --help | head -5

# Run as non-root user
RUN useradd -m echouser
USER echouser

