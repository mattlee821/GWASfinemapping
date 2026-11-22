FROM mambaorg/micromamba:2.0.5

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps binutils file \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/gwasfinemapping /opt/gwasfinemapping-r-lib \
    && chown -R $MAMBA_USER:$MAMBA_USER /opt/gwasfinemapping /opt/gwasfinemapping-r-lib \
    && printf 'export PATH="/opt/conda/bin:${PATH}"\nexport JAVA_HOME="/opt/conda/lib/jvm"\nexport R_LIBS_USER="/opt/gwasfinemapping-r-lib/current"\nexport R_PROFILE_USER=""\n' > /etc/profile.d/gwasfm-conda.sh
USER $MAMBA_USER

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements-tools.txt /tmp/requirements-tools.txt
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y
RUN PYTHONNOUSERSITE=1 micromamba run -n base python -m pip install --upgrade pip \
    && PYTHONNOUSERSITE=1 micromamba run -n base python -m pip install -r /tmp/requirements-tools.txt \
    && PYTHONNOUSERSITE=1 micromamba run -n base python -m pip install --no-deps https://github.com/jordanero/MultiSuSiE/archive/refs/heads/main.zip

WORKDIR /opt/gwasfinemapping
ENV RENV_PROJECT=/opt/gwasfinemapping \
    R_PROFILE_USER=/opt/gwasfinemapping/.Rprofile \
    RENV_PATHS_LIBRARY=/opt/gwasfinemapping-r-lib \
    RENV_CONFIG_PROMPT=FALSE \
    RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

COPY --chown=$MAMBA_USER:$MAMBA_USER renv.lock .Rprofile /opt/gwasfinemapping/
COPY --chown=$MAMBA_USER:$MAMBA_USER renv/ /opt/gwasfinemapping/renv/
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/install_carma.R /opt/gwasfinemapping/scripts/install_carma.R
RUN set -eux; \
    GWASFM_FORCE_RENV=1 Rscript -e 'renv::restore(prompt = FALSE)'; \
    R_LIB="$(GWASFM_FORCE_RENV=1 Rscript -e 'cat(.libPaths()[1])')"; \
    test -n "$R_LIB" && test -d "$R_LIB"; \
    R_LIBS="$R_LIB" R_LIBS_USER="$R_LIB" R_PROFILE_USER="" Rscript scripts/install_carma.R "$R_LIB"; \
    ln -sfn "$R_LIB" /opt/gwasfinemapping-r-lib/current

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/gwasfinemapping
RUN set -eux; \
    R_LIB="$(readlink -f /opt/gwasfinemapping-r-lib/current)"; \
    R_LIBS="$R_LIB" R_LIBS_USER="$R_LIB" R_PROFILE_USER="" R CMD INSTALL --no-multiarch --library="$R_LIB" .; \
    find src \( -name '*.o' -o -name '*.so' \) -delete 2>/dev/null || true; \
    R_LIBS_USER=/opt/gwasfinemapping-r-lib/current R_PROFILE_USER="" Rscript -e 'stopifnot(requireNamespace("susieR", quietly = TRUE), requireNamespace("CARMA", quietly = TRUE), requireNamespace("jsonlite", quietly = TRUE), requireNamespace("testthat", quietly = TRUE), requireNamespace("gwasfinemapping", quietly = TRUE))'

USER $MAMBA_USER
ENV PATH=/opt/conda/bin:$PATH \
    JAVA_HOME=/opt/conda/lib/jvm \
    PYTHONNOUSERSITE=1 \
    GWASFM_IN_CONTAINER=1 \
    RENV_PROJECT=/opt/gwasfinemapping \
    RENV_PATHS_LIBRARY=/opt/gwasfinemapping-r-lib \
    R_LIBS_USER=/opt/gwasfinemapping-r-lib/current \
    R_PROFILE_USER="" \
    RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
