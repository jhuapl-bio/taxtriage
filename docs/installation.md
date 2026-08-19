# Installation

TaxTriage requires two primary dependencies:

1. **Nextflow** — the workflow engine
2. **Docker** or **Singularity** — for containerized module execution

---

## 1. Install Nextflow

Follow instructions at [nf-co.re/docs/usage/installation](https://nf-co.re/docs/usage/installation) or run the following in a WSL2, native Linux, or macOS terminal:

```bash
# Verify Java v11+ is installed
java -version

# Download and install Nextflow
curl -fsSL get.nextflow.io | bash
```

Move the binary to your `$PATH`:

```bash
# Personal path (no sudo required)
mv nextflow ~/bin/

# System-wide (requires sudo)
sudo mv nextflow /usr/local/bin
```

Verify the installation:

```bash
nextflow -v
```

> **HPC users:** Make sure `nextflow` is in your `$PATH` if it's not globally available on your cluster.

---

## 2. Install a Container Runtime

You only need **one** of the following. Docker is recommended for local workstations; Singularity is common on HPCs.

### A. Docker (Recommended)

Follow the OS-specific instructions at [docs.docker.com/engine/install](https://docs.docker.com/engine/install/).

> **Windows users (WSL2):** Install Docker Desktop for Windows — it will be available automatically inside your WSL2 environment.

### B. Singularity (HPC)

Follow the [Singularity installation guide](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

> **Important:** Only **SingularityCE v4+** is supported. Consult your HPC IT team if Singularity is not already available as a module.

When using Singularity, it is strongly recommended to pre-download containers before running the pipeline offline:

```bash
nf-core download https://github.com/jhuapl-bio/taxtriage --singularity-cache-only
```

Set the cache directory to avoid re-downloading images:

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity/cache
```

---

## 3. Optional: Limit Nextflow JVM Memory

In some environments, the Nextflow Java Virtual Machine requests excessive memory. Add the following to your `~/.bashrc` or `~/.bash_profile`:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

---

## Next Steps

- Proceed to [Quick Start](quick-start.md) to run a test pipeline.
- See [Running the Pipeline](running-the-pipeline.md) for your own data.
