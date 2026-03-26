#!/bin/bash

##############################################
# ROOT → PARQUET converter for RGD data
# Style matched to HIPO → ROOT script
##############################################

# --------- CONFIGURATION ---------
BASE_INPUT="/volatile/clas12/suman/00_RGD_Analysis/data/experimental/root"
BASE_OUTPUT="/volatile/clas12/suman/00_RGD_Analysis/data/experimental/parquet"
TEST_OUT_DIR="/lustre24/expphy/volatile/clas12/suman/00_RGD_Analysis/test_runs/01_parquet"
PYTHON_SCRIPT="root_2_parquet.py"
MAX_JOBS=16
# -------------------------------

if [ -z "$BASH_VERSION" ]; then
    echo "Please run this script with bash."
    exit 1
fi

mkdir -p "${BASE_OUTPUT}"
mkdir -p "${TEST_OUT_DIR}"

##############################################
# 0. CHOOSE MODE
##############################################
echo "=========================================================="
echo "            ROOT → PARQUET CONVERSION TOOL"
echo "=========================================================="
echo "Choose Operation Mode:"
echo "  1) Production (Full files, Smart Skipping)"
echo "  2) Validation (50k limit, Always Overwrite)"
read -rp "Choice [1-2]: " MODE_CHOICE

ENTRY_STOP_ARG=""
[[ "$MODE_CHOICE" == "2" ]] && ENTRY_STOP_ARG="--entry-stop 50000"

##############################################
# 1. Choose TARGET
##############################################
echo -e "\nLogical Targets for RGD Analysis:"
echo "  [1] Cu (Copper) | [2] Sn (Tin) | [3] CxC (Carbon) | [4] LD2 (Deuterium)"
echo "----------------------------------------------------------"
read -rp "Enter Target Index: " T_IDX

case $T_IDX in
    1) TARGET="Cu";  FOLDER="CuSn" ;;
    2) TARGET="Sn";  FOLDER="CuSn" ;;
    3) TARGET="CxC"; FOLDER="CxC"  ;;
    4) TARGET="LD2"; FOLDER="LD2"  ;;
    *) echo "ERROR: Invalid selection"; exit 1 ;;
esac

INPUT_TARGET_DIR="${BASE_INPUT}/${FOLDER}"
if [[ ! -d "${INPUT_TARGET_DIR}" ]]; then
    echo "ERROR: Input target directory does not exist: ${INPUT_TARGET_DIR}"
    exit 1
fi

##############################################
# 2. Choose RUN(S)
##############################################
mapfile -t ALL_RUN_DIRS < <(find "${INPUT_TARGET_DIR}" -maxdepth 1 -type d -not -path "${INPUT_TARGET_DIR}" -printf "%f\n" | sort)

echo "------------------------------------------------"
echo "Found ${#ALL_RUN_DIRS[@]} run directories in folder $FOLDER."
echo "Choose Run Selection Mode:"
echo "  1) Single Run (Interactive file selection)"
echo "  2) Multiple Runs by INDEX (e.g., 1 5 10-15)"
echo "  3) ALL Runs"
echo "------------------------------------------------"
read -rp "Choice [1-3]: " RUN_MODE

SELECTED_RUNS=()

if [[ "$RUN_MODE" == "1" || "$RUN_MODE" == "2" ]]; then
    echo "Available Runs:"
    for i in "${!ALL_RUN_DIRS[@]}"; do
        printf "  [%3d] %s" "$((i+1))" "${ALL_RUN_DIRS[$i]}"
        [[ $(( (i+1) % 4 )) -eq 0 ]] && echo ""
    done
    echo -e "\n"

    if [[ "$RUN_MODE" == "1" ]]; then
        read -rp "Enter Index of the run: " R_IDX
        SELECTED_RUNS+=("${ALL_RUN_DIRS[$((R_IDX-1))]}")
    else
        read -rp "Enter indices/ranges (e.g., 1 2 5-10): " R_INPUT
        for item in $R_INPUT; do
            if [[ $item =~ ^([0-9]+)-([0-9]+)$ ]]; then
                for (( j=${BASH_REMATCH[1]}; j<=${BASH_REMATCH[2]}; j++ )); do
                    SELECTED_RUNS+=("${ALL_RUN_DIRS[$((j-1))]}")
                done
            else
                # FIXED: Changed ALL_RUNS to ALL_RUN_DIRS
                SELECTED_RUNS+=("${ALL_RUN_DIRS[$((item-1))]}")
            fi
        done
    fi
elif [[ "$RUN_MODE" == "3" ]]; then
    SELECTED_RUNS=("${ALL_RUN_DIRS[@]}")
fi

##############################################
# 3. Collect Files & Show Summary
##############################################
FILES_TO_PROCESS=()
declare -A FILE_TO_RUN 
TOTAL_COUNT=0

echo -e "\n--- SELECTION SUMMARY ($TARGET) ---"
printf "%-10s | %-15s | %-10s\n" "Index" "Run Number" "File Count"
echo "------------------------------------------------"

for RUN_NUM in "${SELECTED_RUNS[@]}"; do
    RUN_DIR="${INPUT_TARGET_DIR}/${RUN_NUM}"
    mapfile -t ROOT_FILES < <(ls "${RUN_DIR}"/*.root 2>/dev/null | sort || true)
    COUNT=${#ROOT_FILES[@]}
    
    # Matching original script indexing for summary display
    IDX="?"
    for i in "${!ALL_RUN_DIRS[@]}"; do
        [[ "${ALL_RUN_DIRS[$i]}" == "$RUN_NUM" ]] && IDX=$((i+1)) && break
    done

    # FIXED: Changed %3d to %3s to handle the "?" string safely
    printf "[%3s]      | %-15s | %-10d\n" "$IDX" "$RUN_NUM" "$COUNT"

    if [[ "$RUN_MODE" == "1" && $COUNT -gt 0 ]]; then
        echo -e "\nFiles in $RUN_NUM:"
        for i in "${!ROOT_FILES[@]}"; do
            printf "  [%3d] %s\n" "$((i+1))" "$(basename "${ROOT_FILES[$i]}")"
        done
        read -rp "Select: 1) Single 2) First N 3) All 4) Range: " F_MODE
        case $F_MODE in
            1) read -rp "Idx: " I; TEMP_LIST=("${ROOT_FILES[$((I-1))]}") ;;
            2) read -rp "N: " N; TEMP_LIST=("${ROOT_FILES[@]:0:N}") ;;
            3) TEMP_LIST=("${ROOT_FILES[@]}") ;;
            4) read -rp "Range (S-E): " R; [[ $R =~ ^([0-9]+)-([0-9]+)$ ]]; S=${BASH_REMATCH[1]}; E=${BASH_REMATCH[2]}
               TEMP_LIST=("${ROOT_FILES[@]:$((S-1)):$((E-S+1))}") ;;
        esac
    else
        TEMP_LIST=("${ROOT_FILES[@]}")
    fi

    for F in "${TEMP_LIST[@]}"; do
        FILES_TO_PROCESS+=("$F")
        FILE_TO_RUN["$F"]="$RUN_NUM"
        ((TOTAL_COUNT++))
    done
done

echo "------------------------------------------------"
echo "TOTAL FILES QUEUED: $TOTAL_COUNT"
echo "------------------------------------------------"
read -rp "Proceed with ROOT to Parquet conversion? (y/n): " CONFIRM
[[ "$CONFIRM" != "y" ]] && exit 0

##############################################
# 4. Smart Processing Logic
##############################################
process_file() {
    local in_file="$1"
    local run_num="${FILE_TO_RUN[$in_file]}"
    local base_name="$(basename "${in_file}" .root)"
    
    # Mode-based directory logic
    if [[ "$MODE_CHOICE" == "2" ]]; then
        local out_dir="${TEST_OUT_DIR}"
    else
        local out_dir="${BASE_OUTPUT}/${TARGET}/${run_num}"
    fi

    # 1. Skip Logic (Enabled only for Production Mode)
    if [[ "$MODE_CHOICE" == "1" ]]; then
        if ls "${out_dir}"/*"${TARGET}"*"${base_name}"*.parquet >/dev/null 2>&1; then
            echo "[SKIP - EXISTS] Run:${run_num} | File:${base_name} | Location: ${out_dir}"
            return 0
        fi
    fi

    mkdir -p "${out_dir}"

    # Execution (Hiding python logs but keeping shell status)
    python3 "${PYTHON_SCRIPT}" \
        --target "${TARGET}" \
        --root-file "${in_file}" \
        --out-dir "${out_dir}" \
        ${ENTRY_STOP_ARG} > /dev/null 2>&1

    if [[ $? -eq 0 ]]; then
        echo "[DONE] Run:${run_num} | Created Parquet for: ${base_name}"
    else
        echo "[ERROR] Run:${run_num} | Failed: ${base_name} | Check python script"
    fi
}

##############################################
# 5. Execution
##############################################
START_TIME=$(date +%s)
echo "Starting conversion (Parallel: $MAX_JOBS)..."

for f in "${FILES_TO_PROCESS[@]}"; do
    process_file "${f}" &
    while (( $(jobs -r | wc -l) >= MAX_JOBS )); do
        sleep 2
    done
done

wait
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "------------------------------------------------"
echo "Workflow Complete!"
echo "Total Time: $((ELAPSED/3600))h $(((ELAPSED%3600)/60))m $((ELAPSED%60))s"
echo "------------------------------------------------"