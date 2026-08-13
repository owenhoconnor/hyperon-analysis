#!/bin/bash

#SBATCH --job-name=Hyperons_analyzer-NEW
#SBATCH --output=/data/ooconnor/slurm_logs/hyperons/new_logic/filt_hyps/lar_%a.out # Where log files for each job end up
#SBATCH --error=/data/ooconnor/slurm_logs/hyperons/new_logic/filt_hyps/lar_%a.err # Where err files for each job end up
#SBATCH --array=0-999%100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=12:00:00

# Define the relevant paths (CHANGE TO YOURS)
INPUT_DIR="/data/sbnd/hyperons_antinumu"
OUTPUT_DIR="/data/ooconnor/sbnd/hyperons/analyzer_output/new_logic/filtered_hyperons/"
WORK_DIR="$HOME/larsoft_dev/srcs/sbndcode/sbndcode/Hyperons"
SETUP_LOCAL="$HOME/larsoft_dev/localProducts_larsoft_v10_21_02_prof_e26/setup"
CONTAINER="/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest"
APPTAINER_BIN="/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer"

# Grab the list of files
FILES=($INPUT_DIR/firstbatch/prodgenie* $INPUT_DIR/secondbatch/prodgenie* $INPUT_DIR/thirdbatch/prodgenie* $INPUT_DIR/fourhtbatch/prodgenie*) # The name here is the wildcard of the files in the input_dir
NUM_FILES=${#FILES[@]}
NUM_JOBS=1000

# Container stuff
$APPTAINER_BIN exec -B /cvmfs,/data,/home,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf \
    --ipc --pid "$CONTAINER" /bin/bash -c "
        source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
        source $SETUP_LOCAL
        mrbsetenv
        mrbslp

        INPUT_FILES=($INPUT_DIR/firstbatch/prodgenie* $INPUT_DIR/secondbatch/prodgenie* $INPUT_DIR/thirdbatch/prodgenie* $INPUT_DIR/fourhtbatch/prodgenie*)

        # Note: We use 'seq' to handle looking for files
        for (( i=$SLURM_ARRAY_TASK_ID; i<$NUM_FILES; i+=$NUM_JOBS )); do


            CURRENT_FILE=\${INPUT_FILES[\$i]}

            echo \"Processing file index \$i: \$CURRENT_FILE\"

            # Store in temp file
            TMP_DIR=$WORK_DIR/tmp_job_\${SLURM_ARRAY_TASK_ID}_\$i
            mkdir -p \$TMP_DIR
            cd \$TMP_DIR

            # Run lar
            lar -c $WORK_DIR/run_analyzeEvents.fcl \
                -s \$CURRENT_FILE \
                -T analyzer_output_\${i}_2026.root

             # Check for errors, rerun if any jobs failed

            LAR_STATUS=\$?

            if [ \$LAR_STATUS -ne 0 ]; then
                echo \"lar failed for file \$i with exit code \$LAR_STATUS. Rerunning...\"

                rm -f analyzer_output_\${i}_2026.root

                lar -c $WORK_DIR/run_analyzeEvents.fcl \
                -s \$CURRENT_FILE \
                -T analyzer_output_\${i}_2026.root
            fi

            # Move and clean (If neccessary)
            mv analyzer_output_\${i}_2026.root $OUTPUT_DIR/
            cd $WORK_DIR
            rm -rf \$TMP_DIR
        done
    "
