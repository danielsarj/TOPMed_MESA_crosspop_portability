#!/bin/bash

declare -a files
files=("Height" "QRS_duration" "C-reactive_protein" "BMI" "Chronic_kidney" "Smoking" "Coffee" "Diastolic_blood_pressure" "Glomerular_filtration_rate" "End-stage_renal_disease" "Fasting_blood_glucose" "Fasting_blood_insulin" "Hemoglobin_A1c" "HDL" "LDL" "Hypertension" "Mean_corpuscular_hemoglobin" "Platelet_count" "PR_interval" "QT_interval" "Systolic_blood_pressure" "Total_cholesterol" "Triglycerides" "TypeII_diabetes" "WBC_count" "Waist-hip_ratio-50" "Waist-hip_ratio-51" "Waist-hip_ratio-52")

for trait in ${files[@]}
do
    for tissue in PBMC
    do
            for pop in AFA EUR HIS CHN
            do
                    bash /home/daraujo1/scratch/TOPMed_MESA/SPrediXcan/PAGE_SPredixcan.sh ${trait} ${tissue} ${pop}
            done
    done
done

for trait in ${files[@]}
do
    for tissue in Mono Tcell
    do
            for pop in AFA EUR HIS
            do
                    bash /home/daraujo1/scratch/TOPMed_MESA/SPrediXcan/PAGE_SPredixcan.sh ${trait} ${tissue} ${pop}
            done
    done
done

