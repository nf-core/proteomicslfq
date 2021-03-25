FROM nfcore/base:1.13.3
LABEL authors="Julianus Pfeuffer, Lukas Heumos, Leon Bichmann, Timo Sachsenberg, Yasset Perez-Riverol" \
      description="Docker image containing all software requirements for the nf-core/proteomicslfq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-proteomicslfq-1.1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-proteomicslfq-1.1.0dev > nf-core-proteomicslfq-1.1.0dev.yml
