#!/bin/bash

#===============================================================================
# SCRIPT: run_pipeline.sh
# PURPOSE: RNA-seq Analysis Pipeline Runner
#
# DESCRIPTION:
# This script orchestrates the complete RNA-seq analysis pipeline by
# submitting jobs in the correct order with appropriate dependencies.
# It provides a convenient way to run the entire analysis or specific
# steps of the pipeline.
#
# PIPELINE STEPS:
# 1. FastQC - Quality control of raw reads
# 2. Trim Galore - Adapter trimming and quality filtering
# 3. STAR - Read alignment to reference genome
# 4. Quantification - Gene expression quantification
# 5. DESeq2 - Differential expression analysis
#
# USAGE:
# ./run_pipeline.sh [options]
#
# OPTIONS:
#   --all           Run complete pipeline (default)
#   --from-step N   Start from step N (1-5)
#   --to-step N     Stop at step N (1-5)
#   --step N        Run only step N
#   --dry-run       Show what would be submitted without submitting
#   --help          Show this help message
#
# EXAMPLES:
# ./run_pipeline.sh                    # Run complete pipeline
# ./run_pipeline.sh --from-step 3      # Start from alignment
# ./run_pipeline.sh --step 1           # Run only FastQC
# ./run_pipeline.sh --to-step 3        # Run up to alignment
# ./run_pipeline.sh --dry-run          # Preview job submissions
#
# DEPENDENCIES:
# - SLURM workload manager
# - All analysis scripts (1_fastqc.sh through 5_deseq2.sh)
# - Required conda environments and software
#===============================================================================

# Default parameters
FROM_STEP=1
TO_STEP=5
DRY_RUN=false
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_color() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

# Function to show help
show_help() {
    cat << EOF
RNA-seq Analysis Pipeline Runner

USAGE:
    ./run_pipeline.sh [options]

OPTIONS:
    --all           Run complete pipeline (default)
    --from-step N   Start from step N (1-5)
    --to-step N     Stop at step N (1-5)
    --step N        Run only step N
    --dry-run       Show what would be submitted without submitting
    --help          Show this help message

PIPELINE STEPS:
    1. FastQC          - Quality control of raw reads
    2. Trim Galore     - Adapter trimming and quality filtering
    3. STAR Alignment  - Read alignment to reference genome
    4. Quantification  - Gene expression quantification
    5. DESeq2          - Differential expression analysis

EXAMPLES:
    ./run_pipeline.sh                    # Run complete pipeline
    ./run_pipeline.sh --from-step 3      # Start from alignment
    ./run_pipeline.sh --step 1           # Run only FastQC
    ./run_pipeline.sh --to-step 3        # Run up to alignment
    ./run_pipeline.sh --dry-run          # Preview job submissions

NOTES:
    - Jobs are submitted with appropriate dependencies
    - Check SLURM queue with: squeue -u \$USER
    - Cancel jobs with: scancel <job_id>
    - Monitor logs in: ${BASE_DIR}/logs/

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --all)
            FROM_STEP=1
            TO_STEP=5
            shift
            ;;
        --from-step)
            FROM_STEP="$2"
            shift 2
            ;;
        --to-step)
            TO_STEP="$2"
            shift 2
            ;;
        --step)
            FROM_STEP="$2"
            TO_STEP="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            print_color $RED "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Validate step numbers
if [[ $FROM_STEP -lt 1 || $FROM_STEP -gt 5 ]]; then
    print_color $RED "Error: FROM_STEP must be between 1 and 5"
    exit 1
fi
if [[ $TO_STEP -lt 1 || $TO_STEP -gt 5 ]]; then
    print_color $RED "Error: TO_STEP must be between 1 and 5"
    exit 1
fi
if [[ $FROM_STEP -gt $TO_STEP ]]; then
    print_color $RED "Error: FROM_STEP cannot be greater than TO_STEP"
    exit 1
fi

# Change to base directory
cd "$BASE_DIR"

# Create logs directory if it doesn't exist
mkdir -p logs

print_color $BLUE "=========================================="
print_color $BLUE "    RNA-seq Analysis Pipeline Runner"
print_color $BLUE "=========================================="
echo
print_color $YELLOW "Pipeline configuration:"
echo "  Base directory: $BASE_DIR"
echo "  Script directory: $SCRIPT_DIR"
echo "  From step: $FROM_STEP"
echo "  To step: $TO_STEP"
echo "  Dry run: $DRY_RUN"
echo

# Define pipeline steps
declare -A STEPS
STEPS[1]="1_fastqc.sh"
STEPS[2]="2_trim.sh"
STEPS[3]="3_align.sh"
STEPS[4]="4_quantify.sh"
STEPS[5]="5_deseq2.sh"

declare -A STEP_NAMES
STEP_NAMES[1]="FastQC Quality Control"
STEP_NAMES[2]="Adapter Trimming"
STEP_NAMES[3]="STAR Alignment"
STEP_NAMES[4]="Gene Quantification"
STEP_NAMES[5]="Differential Expression"

# Check if all required scripts exist
print_color $YELLOW "Checking required scripts..."
ALL_SCRIPTS_EXIST=true
for i in $(seq $FROM_STEP $TO_STEP); do
    SCRIPT_PATH="$SCRIPT_DIR/${STEPS[$i]}"
    if [[ -f "$SCRIPT_PATH" ]]; then
        print_color $GREEN "  ✓ ${STEPS[$i]} found"
    else
        print_color $RED "  ✗ ${STEPS[$i]} missing"
        ALL_SCRIPTS_EXIST=false
    fi
done

if [[ "$ALL_SCRIPTS_EXIST" == "false" ]]; then
    print_color $RED "Error: Some required scripts are missing. Please check the scripts directory."
    exit 1
fi

echo

# Check required files and directories
print_color $YELLOW "Checking required files and directories..."

# Check sample configuration
if [[ -f "config/samples.txt" ]]; then
    SAMPLE_COUNT=$(wc -l < config/samples.txt)
    print_color $GREEN "  ✓ Sample configuration found ($SAMPLE_COUNT samples)"
else
    print_color $RED "  ✗ config/samples.txt missing"
    exit 1
fi

# Check design file
if [[ -f "config/design.txt" ]]; then
    print_color $GREEN "  ✓ Design file found"
else
    print_color $RED "  ✗ config/design.txt missing"
    exit 1
fi

# Check raw data directory (if starting from step 1 or 2)
if [[ $FROM_STEP -le 2 ]]; then
    RAW_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/download/TES_RNA/TES"
    if [[ -d "$RAW_DIR" ]]; then
        print_color $GREEN "  ✓ Raw data directory found"
    else
        print_color $RED "  ✗ Raw data directory missing: $RAW_DIR"
        exit 1
    fi
fi

echo

# Function to submit job
submit_job() {
    local step=$1
    local script=${STEPS[$step]}
    local step_name=${STEP_NAMES[$step]}
    local dependency=$2

    print_color $YELLOW "Step $step: $step_name" >&2

    if [[ "$DRY_RUN" == "true" ]]; then
        if [[ -n "$dependency" ]]; then
            echo "  [DRY RUN] Would submit: sbatch --dependency=afterok:$dependency scripts/$script" >&2
        else
            echo "  [DRY RUN] Would submit: sbatch scripts/$script" >&2
        fi
        echo "dummy_job_id"
        return 0
    else
        if [[ -n "$dependency" ]]; then
            JOB_ID=$(sbatch --dependency=afterok:$dependency --parsable scripts/$script 2>/dev/null)
        else
            JOB_ID=$(sbatch --parsable scripts/$script 2>/dev/null)
        fi

        if [[ $? -eq 0 && -n "$JOB_ID" ]]; then
            print_color $GREEN "  ✓ Job submitted: $JOB_ID" >&2
            echo "$JOB_ID"
        else
            print_color $RED "  ✗ Failed to submit job" >&2
            return 1
        fi
    fi
}

# Submit jobs with dependencies
print_color $BLUE "Submitting pipeline jobs..."
echo

PREV_JOB_ID=""

for i in $(seq $FROM_STEP $TO_STEP); do
    if [[ $i -eq $FROM_STEP ]]; then
        # First job has no dependency
        PREV_JOB_ID=$(submit_job $i "")
    else
        # Subsequent jobs depend on previous job
        PREV_JOB_ID=$(submit_job $i "$PREV_JOB_ID")
    fi

    if [[ $? -ne 0 ]]; then
        print_color $RED "Pipeline submission failed at step $i"
        exit 1
    fi
done

echo

if [[ "$DRY_RUN" == "true" ]]; then
    print_color $YELLOW "Dry run completed. No jobs were actually submitted."
else
    print_color $GREEN "Pipeline submitted successfully!"
    echo
    print_color $YELLOW "Useful commands:"
    echo "  Check job status:     squeue -u \$USER"
    echo "  Cancel all jobs:      scancel -u \$USER"
    echo "  Monitor logs:         tail -f logs/*.out"
    echo "  View job details:     scontrol show job <job_id>"
fi

echo
print_color $BLUE "Pipeline runner completed."


---

# ● You can run the pipeline starting from step 2 (trimming) using the --from-step option:

# ./scripts/run_pipeline.sh --from-step 2

# This will skip the FastQC step (step 1) and start directly with trimming. Since the trimming appears to have
# completed successfully based on the files we saw in the output directory, you might want to start from step 3
# (alignment):

# ./scripts/run_pipeline.sh --from-step 3

# Or if you want to see what would be submitted first, use a dry run:

# ./scripts/run_pipeline.sh --from-step 3 --dry-run

# The pipeline runner supports these options:
# - --from-step N - Start from step N
# - --to-step N - Stop at step N
# - --step N - Run only step N
# - --dry-run - Preview without submitting

# Since I fixed the file naming issues in the alignment script, starting from step 3 should work correctly with the
# trimmed files that are already generated.