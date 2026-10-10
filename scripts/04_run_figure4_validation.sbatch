#!/bin/bash
#SBATCH --job-name=sparq_fig4
#SBATCH --output=/ourdisk/hpc/wei_r_chen/rahhulse100/dont_archive/SPARQ_figure4_validation/logs/figure4_%j.out
#SBATCH --error=/ourdisk/hpc/wei_r_chen/rahhulse100/dont_archive/SPARQ_figure4_validation/logs/figure4_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4

# Run Figure 4 analyses sequentially, retaining every intermediate output.
set -euo pipefail

ENV_ROOT="/ourdisk/hpc/wei_r_chen/rahhulse100/dont_archive/SPARQ_figure4_envs"
OUT_ROOT="/ourdisk/hpc/wei_r_chen/rahhulse100/dont_archive/SPARQ_figure4_validation"
DATA_ROOT="/ourdisk/hpc/wei_r_chen/rahhulse100/dont_archive/from_scratch_rahhulse100/spatialPDAC_copy/SPARQ/validation_datasets"
R_ENV="${ENV_ROOT}/fig4_r_bioc"
PY_ENV="${ENV_ROOT}/fig4_py"
SPAGCN_ENV="${ENV_ROOT}/fig4_spagcn"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUN_ID="run_${SLURM_JOB_ID}"
RUN_DIR="${OUT_ROOT}/${RUN_ID}"

mkdir -p "${RUN_DIR}"/{inputs,raw,standardized,sparq,logs,manifests,software,scripts}
cp "${SCRIPT_DIR}"/{01_prepare_figure4_inputs.py,02_run_python_methods.py,03_run_r_methods.R,"$(basename "$0")"} "${RUN_DIR}/scripts/"

if [ -f "${HOME}/.bashrc" ]; then source "${HOME}/.bashrc" || true; fi
source "$(conda info --base)/etc/profile.d/conda.sh"

for env in "${R_ENV}" "${PY_ENV}" "${SPAGCN_ENV}"; do
  [ -d "${env}" ] || { echo "Missing environment: ${env}" >&2; exit 2; }
done

{
  echo -e "key\tvalue"
  echo -e "run_id\t${RUN_ID}"
  echo -e "started\t$(date -Is)"
  echo -e "data_root\t${DATA_ROOT}"
  echo -e "r_environment\t${R_ENV}"
  echo -e "python_environment\t${PY_ENV}"
  echo -e "spagcn_environment\t${SPAGCN_ENV}"
  echo -e "slurm_job_id\t${SLURM_JOB_ID}"
  echo -e "cpus\t${SLURM_CPUS_PER_TASK}"
} > "${RUN_DIR}/manifests/run_manifest.tsv"

git clone --depth 1 https://github.com/rahhulse100/SPARQ1.git "${RUN_DIR}/software/SPARQ1" \
  > "${RUN_DIR}/logs/00_clone_sparq.log" 2>&1
git -C "${RUN_DIR}/software/SPARQ1" rev-parse HEAD > "${RUN_DIR}/manifests/sparq_git_commit.txt"
mkdir -p "${RUN_DIR}/software/R_library"
conda run --no-capture-output -p "${R_ENV}" \
  R CMD INSTALL -l "${RUN_DIR}/software/R_library" "${RUN_DIR}/software/SPARQ1" \
  > "${RUN_DIR}/logs/00_install_sparq.log" 2>&1

conda run --no-capture-output -p "${PY_ENV}" python "${SCRIPT_DIR}/01_prepare_figure4_inputs.py" \
  --data-root "${DATA_ROOT}" --output-root "${RUN_DIR}" --expected-sections 14 \
  > "${RUN_DIR}/logs/01_prepare_inputs.log" 2>&1

conda run --no-capture-output -p "${PY_ENV}" python "${SCRIPT_DIR}/02_run_python_methods.py" \
  --run-dir "${RUN_DIR}" --python-env "${PY_ENV}" --spagcn-env "${SPAGCN_ENV}" --iterations 50 --retention 0.75 --seed 20261010 \
  > "${RUN_DIR}/logs/02_python_methods.log" 2>&1

conda run --no-capture-output -p "${R_ENV}" Rscript --vanilla "${SCRIPT_DIR}/03_run_r_methods.R" \
  "${RUN_DIR}" "${R_ENV}" 50 0.75 20261010 \
  > "${RUN_DIR}/logs/03_r_methods_and_sparq.log" 2>&1

date -Is > "${RUN_DIR}/manifests/completed_at.txt"
echo "Figure 4 workflow completed: ${RUN_DIR}"
