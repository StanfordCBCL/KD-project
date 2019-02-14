#!/bin/bash 

# python prelim_analysis.py --suff p1 --vtk --vtu --source AneurysmGeneration/models/SKD0050/right_asi_2/SKD0050_baseline_model_modified_p1.vtp --results Artificial/RCA/ASI2/p1/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/RCA/ASI2/ --suff p1 --shape ASI2
# python prelim_analysis.py --suff p5 --vtk --vtu --source AneurysmGeneration/models/SKD0050/right_asi_2/SKD0050_baseline_model_modified_p5.vtp --results Artificial/RCA/ASI2/p5/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/RCA/ASI2/ --suff p5 --shape ASI2

# python prelim_analysis.py --suff p1 --vtk --vtu --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p1.vtp --results Artificial/RCA/ASI4/p1/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/RCA/ASI4/ --suff p1 --shape ASI4
# python prelim_analysis.py --suff p5 --vtk --vtu --source AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p5.vtp --results Artificial/RCA/ASI4/p5/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/RCA/ASI4/ --suff p5 --shape ASI4

python prelim_analysis.py --suff p1 --vtk --vtu --source AneurysmGeneration/models/SKD0050/right_asi_6/SKD0050_baseline_model_modified_p1.vtp --results Artificial/RCA/ASI6/p1/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/RCA/ASI6/ --suff p1 --shape ASI6
python prelim_analysis.py --suff p5 --vtk --vtu --source AneurysmGeneration/models/SKD0050/right_asi_6/SKD0050_baseline_model_modified_p5.vtp --results Artificial/RCA/ASI6/p5/all_results.vtp --mapping --post --clip  --outdir clipped_results_short/RCA/ASI6/ --suff p5 --shape ASI6
