# Use an official Python runtime as a parent image
FROM python:3.10-slim

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# Install Python dependencies
RUN pip install --no-cache-dir \
    dendropy==4.6.4 \
    toytree==3.0.4 \
    toyplot==1.0.3 \
    pandas \
    numpy \
    seaborn==0.13.2 \
    matplotlib \
    biopython==1.79 \
    scikit-learn==1.5.1
