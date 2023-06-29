# Use the lightweight Alpine image as the base
FROM alpine:latest

# Add Python and pip to the image
RUN apk add --no-cache python3 py3-pip gzip gcc

# Ensure Python & pip are up to date
RUN python3 -m ensurepip && \
    pip3 install --upgrade pip setuptools

# Install the required Python modules
RUN pip3 install biopython

# Install BBMap
RUN apk add --no-cache wget bash
RUN mkdir /bbmap && cd /bbmap
RUN wget https://sourceforge.net/projects/bbmap/files/latest/download -O bbmap.tar.gz && \
    tar -xzf bbmap.tar.gz && \
    rm bbmap.tar.gz

# Add bbmap to path
ENV PATH="/bbmap:${PATH}"

# Install PyPy
RUN apk add --no-cache pypy3
RUN pypy3 -m ensurepip
RUN pypy3 -m pip install biopython argparse

# Default command to execute when the container starts
CMD [ "bash" ]
