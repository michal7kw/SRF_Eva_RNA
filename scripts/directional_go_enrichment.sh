#!/bin/bash
#SBATCH --job-name=dir_go
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/directional_go.out
#SBATCH --error=logs/directional_go.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "DIRECTIONAL GO ENRICHMENT (UP vs DOWN)"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA

# Create logs directory if missing
mkdir -p logs

echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

echo "Running directional GO enrichment..."
Rscript scripts/directional_go_enrichment.R

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Analysis completed successfully!"
    echo "=========================================="
    echo "Results: results/06_directional_go/"
else
    echo "ERROR: Analysis failed!"
    exit 1
fi
