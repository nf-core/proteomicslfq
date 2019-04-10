FROM nfcore/base
LABEL authors="The Heumos Brothers - Simon and Lukas" \
      description="Docker image containing all requirements for nf-core/proteomicslfq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-proteomicslfq-1.0dev/bin:$PATH
