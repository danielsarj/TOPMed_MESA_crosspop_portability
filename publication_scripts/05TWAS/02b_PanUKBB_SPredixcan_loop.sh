#!/bin/bash

declare -a files

files=("BMI_calcuated" "BMI_estimated" "Coffee" "C-reactive_protein" "Diastolic_blood_pressure_auto" "Diastolic_blood_pressure_manual" "End-stage_renal_disease" "Fasting_blood_glucose_impaired" "Fasting_blood_glucose" "Glomerular_filtration_rate_cystain_C" "Glomerular_filtration_rate_serum_creatinine_and_cystain_C" "Glomerular_filtration_rate_serum_creatinine" "HDL" "Height_sitting" "Height_standing" "Hemoglobin_A1c" "Hypertension_noncancer" "Hypertension" "LDL" "Mean_corpuscular_hemoglobin" "Platelet" "PR_interval" "QRS_duration" "Smoking" "Systolic_blood_pressure_auto" "Systolic_blood_pressure_manual" "Total_cholesterol" "Triglycerides" "TypeII_diabetes" "Waist-hip_ratio_hip_circumference" "Waist-hip_ratio_waist_circumference" "WBC")

for trait in ${files[@]}
do
    for tissue in PBMC
    do
	for pop in AFR EUR NAM EAS
	do
		for summary_pop in meta
		do
                    bash /home/daraujo1/scratch/TOPMed_MESA/SPrediXcan/PanUKBB_SPredixcan.sh ${trait} ${tissue} ${pop} ${summary_pop}
            	done
	done
    done
done

for trait in ${files[@]}
do
    for tissue in Mono Tcell
    do
	for pop in AFR EUR NAM
	do
		for summary_pop in meta
		do
                    bash /home/daraujo1/scratch/TOPMed_MESA/SPrediXcan/PanUKBB_SPredixcan.sh ${trait} ${tissue} ${pop} ${summary_pop}
            	done
	done
    done
done

for trait in ${files[@]}
do
    for tissue in PBMC Mono Tcell
    do
	for pop in AFR
	do
		for summary_pop in AFR
		do
                    bash /home/daraujo1/scratch/TOPMed_MESA/SPrediXcan/PanUKBB_SPredixcan.sh ${trait} ${tissue} ${pop} ${summary_pop}
            	done
	done
    done
done

