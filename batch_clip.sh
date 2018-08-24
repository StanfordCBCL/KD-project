#!/bin/bash 

# proximal result clipping
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p1.vtp  --results Artificial/RCA/p1/all_results.vtp --mapping --post --suff p1 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p2.vtp  --results Artificial/RCA/p2/all_results.vtp --mapping --post --suff p2 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p3.vtp  --results Artificial/RCA/p3/all_results.vtp --mapping --post --suff p3 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p4.vtp  --results Artificial/RCA/p4/all_results.vtp --mapping --post --suff p4 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p5.vtp  --results Artificial/RCA/p5/all_results.vtp --mapping --post --suff p5 --clip

# # medial results
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_m1.vtp  --results Artificial/RCA/m1/all_results.vtp --mapping --post --suff m1 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_m2.vtp  --results Artificial/RCA/m2/all_results.vtp --mapping --post --suff m2 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_m3.vtp  --results Artificial/RCA/m3/all_results.vtp --mapping --post --suff m3 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_m4.vtp  --results Artificial/RCA/m4/all_results.vtp --mapping --post --suff m4 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_m5.vtp  --results Artificial/RCA/m5/all_results.vtp --mapping --post --suff m5 --clip

# # distal results
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_d1.vtp  --results Artificial/RCA/d1/all_results.vtp --mapping --post --suff d1 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_d2.vtp  --results Artificial/RCA/d2/all_results.vtp --mapping --post --suff d2 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_d3.vtp  --results Artificial/RCA/d3/all_results.vtp --mapping --post --suff d3 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_d4.vtp  --results Artificial/RCA/d4/all_results.vtp --mapping --post --suff d4 --clip
# python prelim_analysis.py --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_d5.vtp  --results Artificial/RCA/d5/all_results.vtp --mapping --post --suff d5 --clip

# right ASI 6
# python prelim_analysis.py --suff p1 --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p1.vtp  --results Artificial/RCA/ASI6/p1/all_results.vtp --mapping --post --clip --outdir clipped_results_short/RCA/ASI6/ --targ_fname /AneurysmGeneration/targets.txt
# python prelim_analysis.py --suff p2 --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p2.vtp  --results Artificial/RCA/ASI6/p2/all_results.vtp --mapping --post --clip --outdir clipped_results_short/RCA/ASI6/ --targ_fname /AneurysmGeneration/targets.txt 
# python prelim_analysis.py --suff p3 --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p3.vtp  --results Artificial/RCA/ASI6/p3/all_results.vtp --mapping --post --clip --outdir clipped_results_short/RCA/ASI6/ --targ_fname /AneurysmGeneration/targets.txt
# python prelim_analysis.py --suff p4 --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p4.vtp  --results Artificial/RCA/ASI6/p4/all_results.vtp --mapping --post --clip --outdir clipped_results_short/RCA/ASI6/ --targ_fname /AneurysmGeneration/targets.txt
# python prelim_analysis.py --suff p5 --vtk --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p5.vtp  --results Artificial/RCA/ASI6/p5/all_results.vtp --mapping --post --clip --outdir clipped_results_short/RCA/ASI6/ --targ_fname /AneurysmGeneration/targets.txt

# LEFT asi 2
# python prelim_analysis.py --suff lad1 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_2/SKD0050_baseline_model_modified_lad1.vtp --results Artificial/LAD/ASI2/lad1/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI2/ --targ_fname /AneurysmGeneration/left_targets.txt
# python prelim_analysis.py --suff lad2 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_2/SKD0050_baseline_model_modified_lad2.vtp --results Artificial/LAD/ASI2/lad2/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI2/ --targ_fname /AneurysmGeneration/left_targets.txt
# python prelim_analysis.py --suff lad3 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_2/SKD0050_baseline_model_modified_lad3.vtp --results Artificial/LAD/ASI2/lad3/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI2/ --targ_fname /AneurysmGeneration/left_targets.txt
# python prelim_analysis.py --suff lad4 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_2/SKD0050_baseline_model_modified_lad4.vtp --results Artificial/LAD/ASI2/lad4/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI2/ --targ_fname /AneurysmGeneration/left_targets.txt
# python prelim_analysis.py --suff lad5 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_2/SKD0050_baseline_model_modified_lad5.vtp --results Artificial/LAD/ASI2/lad5/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI2/ --targ_fname /AneurysmGeneration/left_targets.txt

# LEFT asi 6 
# python prelim_analysis.py --suff lad1 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_6/SKD0050_baseline_model_modified_lad1.vtp --results Artificial/LAD/ASI6/lad1/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI6/ --targ_fname /AneurysmGeneration/left_targets.txt
python prelim_analysis.py --suff lad2 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_6/SKD0050_baseline_model_modified_lad2.vtp --results Artificial/LAD/ASI6/lad2/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI6/ --targ_fname /AneurysmGeneration/left_targets.txt
python prelim_analysis.py --suff lad3 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_6/SKD0050_baseline_model_modified_lad3.vtp --results Artificial/LAD/ASI6/lad3/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI6/ --targ_fname /AneurysmGeneration/left_targets.txt
python prelim_analysis.py --suff lad4 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_6/SKD0050_baseline_model_modified_lad4.vtp --results Artificial/LAD/ASI6/lad4/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI6/ --targ_fname /AneurysmGeneration/left_targets.txt
python prelim_analysis.py --suff lad5 --vtk --source AneurysmGeneration/models/SKD0050/left_asi_6/SKD0050_baseline_model_modified_lad5.vtp --results Artificial/LAD/ASI6/lad5/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/LAD/ASI6/ --targ_fname /AneurysmGeneration/left_targets.txt

