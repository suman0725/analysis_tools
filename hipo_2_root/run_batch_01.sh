onfiguration ---
INPUT_DIR="/volatile/clas12/ouillon/RGD/rgd_simulation_pion_7/recon"
OUTPUT_DIR="/w/hallb-scshelf2102/clas12/suman/c12/ld2_mat_1M"
MACRO_NAME="alert_to_root.cxx"
LOG_FILE="processed_files.txt"
START_INDEX=0
END_INDEX=19
MAX_JOBS=20 # Process up to 20 files simultaneously

# Clear a temporary log file for this run
TEMP_LOG="${LOG_FILE}.temp"
> "${TEMP_LOG}"

# --- Function to run the conversion for a single file ---
process_file() {
    INDEX=$1
    INPUT_FILE="rec_events_lD2_${INDEX}.hipo"
    OUTPUT_FILE="mat_rec_ld2_${INDEX}.root"
    INPUT_PATH="${INPUT_DIR}/${INPUT_FILE}"
    OUTPUT_PATH="${OUTPUT_DIR}/${OUTPUT_FILE}"

    echo "Starting conversion for ${INPUT_FILE}..."

    # The ROOT command to execute your analysis macro
    root -b -q "${MACRO_NAME}+(\"${INPUT_PATH}\",\"${OUTPUT_PATH}\")"

    # Check exit status of the ROOT command
    if [ $? -eq 0 ]; then
        echo "${INPUT_FILE}" >> "${TEMP_LOG}"
        echo "Finished ${INPUT_FILE} successfully."
    else
        echo "ERROR: Failed to process ${INPUT_FILE}. Check ROOT output above."
    fi
}

# --- Main loop to launch jobs ---
echo "Launching jobs for files ${START_INDEX} through ${END_INDEX}..."

for i in $(seq ${START_INDEX} ${END_INDEX}); do
    process_file ${i} &
    
    # Check if we have hit the parallel job limit and wait if necessary
    while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do
        sleep 5
    done
done

# Wait for all background jobs to finish
wait

# --- Final Logging and Cleanup ---
echo "All conversion jobs completed."

# Sort the temporary log numerically and append to the final log file
if [ -f "${TEMP_LOG}" ]; then
    sort -V "${TEMP_LOG}" >> "${LOG_FILE}"
    rm "${TEMP_LOG}"
    echo "Successfully processed files listed in ${LOG_FILE}."
else
    echo "No files were successfully processed."
fi
