# Dockerfile for OPM Flow testing
#
# Builds opm-simulators with all dependencies and test data so that an AI agent
# (or any user) can run the full test suite with a single command.
#
# Build:
#   docker build -t opm-flow-tests .
#
# Run all tests (parallel, with timing):
#   docker run opm-flow-tests ctest --parallel $(nproc) --output-on-failure --test-output-size-passed 65536
#
# Run tests and show per-test timing:
#   docker run opm-flow-tests ctest --parallel $(nproc) --output-on-failure --verbose
#
# List available tests:
#   docker run opm-flow-tests ctest --show-only
#
# Run a specific test:
#   docker run opm-flow-tests ctest -R <test_name> --output-on-failure
#
# Run only unit tests (fast, ~2-5 min):
#   docker run opm-flow-tests ctest -L unit --parallel $(nproc) --output-on-failure
#
# Run only regression tests (~30-60 min):
#   docker run opm-flow-tests ctest -L regression --parallel $(nproc) --output-on-failure
#
# Interactive shell (for debugging or exploration):
#   docker run -it opm-flow-tests bash
#
FROM ubuntu:24.04 AS base

ARG DEBIAN_FRONTEND=noninteractive
ARG OPM_VERSION=master
ARG BUILD_JOBS=4

# ── System dependencies ──────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git ca-certificates ccache \
    libboost-all-dev \
    libdune-common-dev libdune-geometry-dev libdune-grid-dev \
    libdune-istl-dev libdune-localfunctions-dev libdune-uggrid-dev \
    libblas-dev liblapack-dev libsuitesparse-dev \
    libopenmpi-dev openmpi-bin \
    pkg-config zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# ── ccache configuration ─────────────────────────────────────────────
ENV CCACHE_DIR=/ccache
ENV CMAKE_C_COMPILER_LAUNCHER=ccache
ENV CMAKE_CXX_COMPILER_LAUNCHER=ccache

WORKDIR /opm

# ── Build opm-common ─────────────────────────────────────────────────
RUN git clone --depth 1 --branch ${OPM_VERSION} \
        https://github.com/OPM/opm-common.git \
    && cmake -S opm-common -B opm-common/build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_TESTING=OFF \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
    && cmake --build opm-common/build -j ${BUILD_JOBS} \
    && cmake --install opm-common/build

# ── Build opm-grid ───────────────────────────────────────────────────
RUN git clone --depth 1 --branch ${OPM_VERSION} \
        https://github.com/OPM/opm-grid.git \
    && cmake -S opm-grid -B opm-grid/build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_TESTING=OFF \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
    && cmake --build opm-grid/build -j ${BUILD_JOBS} \
    && cmake --install opm-grid/build

# ── Clone opm-tests (reference test data) ────────────────────────────
RUN git clone --depth 1 --branch ${OPM_VERSION} \
        https://github.com/OPM/opm-tests.git

# ── Build opm-simulators with tests ─────────────────────────────────
COPY . /opm/opm-simulators

# Set MPIEXEC_MAX_NUMPROCS high so all parallel tests (including 8-proc
# ones) are registered even when building on a machine with fewer cores.
RUN cmake -S opm-simulators -B opm-simulators/build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DOPM_TESTS_ROOT=/opm/opm-tests \
        -DMPIEXEC_MAX_NUMPROCS=16 \
    && cmake --build opm-simulators/build -j ${BUILD_JOBS}

# ── Runtime configuration ────────────────────────────────────────────
WORKDIR /opm/opm-simulators/build
ENV OPM_TESTS_ROOT=/opm/opm-tests

# Allow running MPI tests as root inside the container
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Allow MPI oversubscription so tests requesting more MPI ranks than
# physical cores (e.g. 8-proc tests on a 2-core runner) still execute.
ENV OMPI_MCA_rmaps_base_oversubscribe=1

# Default: list available tests
CMD ["ctest", "--show-only"]
