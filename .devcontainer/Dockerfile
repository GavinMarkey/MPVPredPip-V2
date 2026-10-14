# Production-oriented image for the Snakemake peptide prediction pipeline.
# Build from the repository root:
#   docker build -f .devcontainer/Dockerfile -t mpvpredpip:prod .

FROM condaforge/miniforge3:24.9.2-0 AS builder

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    CONDA_ENV=/opt/conda/envs/mpvpredpip-v2

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        gawk \
        git \
        tcsh \
        unzip \
        wget \
    && rm -rf /var/lib/apt/lists/*

COPY workflow/envs/core.yaml /tmp/core.yaml
RUN conda env create -f /tmp/core.yaml \
    && conda run -n mpvpredpip-v2 python -m pip install --no-cache-dir \
        "pytest>=7.4" "ruff>=0.4" \
    && conda clean --all --yes \
    && rm -f /tmp/core.yaml

FROM condaforge/miniforge3:24.9.2-0 AS runtime

ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    CONDA_DEFAULT_ENV=mpvpredpip-v2 \
    PATH=/opt/conda/envs/mpvpredpip-v2/bin:/opt/conda/bin:${PATH} \
    PYTHONPATH=/workspaces/MPVPredPip/src \
    SNAKEMAKE_CONDA_PREFIX=/opt/snakemake-conda

# curl is not optional here: setup/install_iedb.sh downloads all five IEDB
# standalones and, before this, had no wget fallback - so a runtime stage
# without curl could not install the predictors at all.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        gawk \
        git \
        tcsh \
        unzip \
        wget \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --gid "${USER_GID}" "${USERNAME}" \
    && useradd --uid "${USER_UID}" --gid "${USER_GID}" --create-home --shell /bin/bash "${USERNAME}" \
    && mkdir -p /opt/snakemake-conda /workspaces/MPVPredPip \
    && chown -R "${USER_UID}:${USER_GID}" /opt/conda /opt/snakemake-conda /workspaces

COPY --from=builder --chown=${USER_UID}:${USER_GID} /opt/conda/envs/mpvpredpip-v2 /opt/conda/envs/mpvpredpip-v2
COPY --chown=${USER_UID}:${USER_GID} src/ /workspaces/MPVPredPip/src/
COPY --chown=${USER_UID}:${USER_GID} workflow/ /workspaces/MPVPredPip/workflow/
COPY --chown=${USER_UID}:${USER_GID} config.yaml pyproject.toml /workspaces/MPVPredPip/

RUN ln -sf /opt/conda/envs/mpvpredpip-v2/bin/python /usr/local/bin/python \
    && ln -sf /opt/conda/envs/mpvpredpip-v2/bin/python /usr/local/bin/python3 \
    && printf 'source /opt/conda/etc/profile.d/conda.sh\nconda activate mpvpredpip-v2\n' > /home/${USERNAME}/.bashrc \
    && chown ${USER_UID}:${USER_GID} /home/${USERNAME}/.bashrc

USER ${USERNAME}
WORKDIR /workspaces/MPVPredPip

CMD ["snakemake", "--use-conda", "--cores", "1"]
