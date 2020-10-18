from openms/executables:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH /opt/conda/bin:$PATH

RUN apt-get install -y wget
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN bash -c "conda init"

# Install the conda environment
COPY environment-dev.yml /
RUN bash -c "conda env create -f /environment-dev.yml && conda clean -a"

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-proteomicslfq-1.0dev/bin:$PATH
CMD [/bin/bash]
