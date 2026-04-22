#!/bin/bash

set -e  # Exit on error



# === Ask for model number ===
read -p "Enter the model number: " model_num

ddir="TTI_3Axes"
prefix="SyntheticCst"

# === Define paths ===
mode_file="/scratch/gpfs/TROMP/we3822/fairhead/Models/${model_num}/ModeList.txt"
mode_dir="/scratch/gpfs/TROMP/we3822/fairhead/Data/JuliaModes/${ddir}/"
output_dir="/scratch/gpfs/TROMP/we3822/fairhead/Models/${model_num}"
maxs_file="${output_dir}/mode_max_s"
params_f90="/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/src/params.f90"
optvani_f90="/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/src/fairhead_optvani.f90"
vani_cu="/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/src/vani.cu"  # ← Update this path

# === Check inputs ===
if [[ ! -f "$mode_file" ]]; then
    echo "❌ Mode file not found: $mode_file"
    exit 1
fi

# === Step 1: Read mode names ===
mapfile -t modes < "$mode_file"
nmodes=${#modes[@]}
echo "✅ Found $nmodes modes."

# === Step 2: Update nmodes in params.f90 ===
sed -i.bak -E "s/(integer,\s*parameter\s*::\s*nmodes\s*=).*/\1 ${nmodes}/" "$params_f90"
echo "✅ Updated nmodes in $params_f90"

# === Step 3: Extract modeNs, modeLs, and dataSmax ===
Ns=()
Ls=()
Smax=()

> "$maxs_file"  # Clear the output file

for mode in "${modes[@]}"; do
    # Extract N (before 'S') and L (after 'S')
    N="${mode%%[ST]*}"
    L="${mode##*[ST]}"
    Ns+=("$N")
    Ls+=("$L")

    # Read line 3 from file
    obs_file="${mode_dir}/${prefix}.${mode}"
    if [[ ! -f "$obs_file" ]]; then
        echo "⚠️  Warning: $obs_file not found. Using -1 for smax."
        smax="-1"
    else
        smax=$(sed -n '3p' "$obs_file" | tr -d '[:space:]')
        if [[ -z "$smax" ]]; then
            echo "⚠️  Warning: Line 3 of $obs_file is empty. Using -1."
            smax="-1"
        fi
    fi
    echo "$smax" >> "$maxs_file"
    Smax+=("$smax")
done

echo "✅ Wrote dataSmax values to $maxs_file"

# === Step 4: Format arrays as comma-separated lists ===
Ns_string=$(IFS=, ; echo "${Ns[*]}")
Ls_string=$(IFS=, ; echo "${Ls[*]}")
Smax_string=$(IFS=, ; echo "${Smax[*]}")

# === Step 5: Update fairhead_optvani.f90 arrays ===
sed -i.bak -E "/:: modeNs\s*=/c\    integer, dimension(nmodes), parameter :: modeNs=(/ $Ns_string /)" "$optvani_f90"
grep -q "modeNs" "$optvani_f90" && echo "✅ modeNs updated" || echo "❌ modeNs sed failed to match"

sed -i -E "/:: modeLs\s*=/c\    integer, dimension(nmodes), parameter :: modeLs=(/ $Ls_string /)" "$optvani_f90"
grep -q "modeLs" "$optvani_f90" && echo "✅ modeLs updated" || echo "❌ modeLs sed failed to match"

sed -i -E "/:: dataSmax\s*=/c\    integer, dimension(nmodes), parameter :: dataSmax=(/ $Smax_string /)" "$optvani_f90"
grep -q "dataSmax" "$optvani_f90" && echo "✅ dataSmax updated" || echo "❌ dataSmax sed failed to match"
echo "✅ Updated modeNs, modeLs, and dataSmax in $optvani_f90"

# === Step 6: Update Lvals in vani.cu ===
if [[ ! -f "$vani_cu" ]]; then
    echo "❌ vani.cu not found: $vani_cu"
    exit 1
fi

# Replace both the array size and values on the active Lvals line
# Matches both commented-out and active lines safely by targeting the uncommented one
sed -i.bak -E "s|^(__constant__ int Lvals\[)[0-9]+(\] = \{)[^}]+(.*)|\\1${nmodes}\\2${Ls_string}\\3|" "$vani_cu"

echo "✅ Updated Lvals[${nmodes}] in $vani_cu"
echo "🎉 All updates completed successfully for model ${model_num}"