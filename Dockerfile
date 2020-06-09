FROM nfcore/base:1.9
LABEL authors="Julianus Pfeuffer, Lukas Heumos, Leon Bichmann, Timo Sachsenberg, Yasset Perez-Riverol" \
      description="Docker image containing all software requirements for the nf-core/proteomicslfq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-proteomicslfq-1.0dev/bin:$PATH

# OpenMS Adapters need the raw jars of Java-based bioconda tools in the PATH. Not the wrappers that conda creates.
RUN cp $(find /opt/conda/envs/nf-core-proteomicslfq-*/share/msgf_plus-*/MSGFPlus.jar -maxdepth 0) $(find /opt/conda/envs/nf-core-proteomicslfq-*/bin/ -maxdepth 0)
RUN cp $(find /opt/conda/envs/nf-core-proteomicslfq-*/share/luciphor2-*/luciphor2.jar -maxdepth 0) $(find /opt/conda/envs/nf-core-proteomicslfq-*/bin/ -maxdepth 0)

# ------------- Parts for dev-only (to have nightly versions of some tools) -------------#

RUN Rscript -e 'if (library("devtools"); install_github("cbielow/PTXQC", build_vignettes = FALSE, dependencies = TRUE)'

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-proteomicslfq-1.0dev > nf-core-proteomicslfq-1.0dev.yml
