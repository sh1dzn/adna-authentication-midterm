# Container for the ancient DNA authentication pipeline.
# NOTE: written but not built during development (no Docker daemon was available);
# the locked conda environment in envs/environment.lock.txt was the tested route.
FROM mambaorg/micromamba:1.5-jammy
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER . /work
CMD ["bash", "scripts/run_test.sh"]
