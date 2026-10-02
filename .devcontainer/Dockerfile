# Dev container image for the peptide-based vaccine candidate prediction pipeline.
#
# The image carries only the orchestrator (Snakemake + the vaxpipe library) and
# the system-level things the predictors need. The heavy, mutually incompatible
# scientific stacks are NOT baked in: Snakemake resolves one conda environment
# per rule from workflow/envs/*.yaml on first run, and those environments are
# cached in a named volume (see docker-compose.yml) so a rebuild of this image
# does not re-resolve them.

FROM condaforge/miniforge3:24.9.2-0

# Non-interactive apt, and unbuffered Python so log output appears in order when
# Snakemake pipes a rule's stdout.
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# System dependencies.
#   git, curl, wget, unzip  - setup/install_iapred.sh and the manual predictor
#                             downloads described in tools/README.md
#   tcsh, gawk              - several IEDB standalone wrappers shell out to them
#   build-essential         - source builds when a conda package has no wheel
#   procps                  - VS Code's server probes the process table
#   less, nano, ca-certificates - ordinary interactive use / TLS
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        gawk \
        git \
        less \
        nano \
        procps \
        tcsh \
        unzip \
        wget \
    && rm -rf /var/lib/apt/lists/*

# A non-root user, so files created in the bind-mounted workspace are not owned
# by root on the host. 1000:1000 matches the default WSL2 / Linux desktop user.
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000
RUN groupadd --gid "${USER_GID}" "${USERNAME}" \
    && useradd --uid "${USER_UID}" --gid "${USER_GID}" --create-home --shell /bin/bash "${USERNAME}" \
    && chown -R "${USER_UID}:${USER_GID}" /opt/conda

# Mount points for the two named volumes (see docker-compose.yml). Docker takes
# a mount point's ownership from the image, so these must exist and be owned by
# the user before the volumes are attached - otherwise they arrive root-owned
# and nothing can write to them.
#
# The conda environment cache deliberately lives OUTSIDE the workspace. Mounting
# a volume at .snakemake/conda would make Docker create .snakemake/ inside the
# bind mount as root, and Snakemake would then fail to create its sibling
# shadow/ and log/ directories.
RUN mkdir -p /opt/snakemake-conda \
    && chown -R "${USER_UID}:${USER_GID}" /opt/snakemake-conda \
    && mkdir -p /workspaces/MPVPredPip \
    && chown -R "${USER_UID}:${USER_GID}" /workspaces

USER ${USERNAME}

# The orchestrator environment. Built from the same file the local-install path
# in README.md uses, so the two stay in step.
COPY --chown=${USER_UID}:${USER_GID} workflow/envs/core.yaml /tmp/core.yaml
RUN conda env create -f /tmp/core.yaml \
    && conda clean --all --yes \
    && rm /tmp/core.yaml

# Development extras. Kept out of core.yaml because that file also describes the
# environment a bare conda install uses to *run* the pipeline, where pytest and
# ruff are not wanted.
RUN conda run -n mpvpredpip-v2 python -m pip install --no-cache-dir \
        "pytest>=7.4" \
        "ruff>=0.4"

# Activate mpvpredpip-v2 for every interactive shell and for `docker exec`.
RUN conda init bash \
    && echo "conda activate mpvpredpip-v2" >> /home/${USERNAME}/.bashrc
ENV PATH=/opt/conda/envs/mpvpredpip-v2/bin:${PATH} \
    CONDA_DEFAULT_ENV=mpvpredpip-v2
# Expose the interpreter at a conventional location as well. VS Code's Python
# extension does not reliably discover environments under /opt/conda/envs, and
# when it finds nothing it offers to install a *second* Python through uv -
# which would have none of this project's dependencies.
USER root
RUN ln -sf /opt/conda/envs/mpvpredpip-v2/bin/python /usr/local/bin/python \
    && ln -sf /opt/conda/envs/mpvpredpip-v2/bin/python /usr/local/bin/python3
USER ${USERNAME}

WORKDIR /workspaces/MPVPredPip
