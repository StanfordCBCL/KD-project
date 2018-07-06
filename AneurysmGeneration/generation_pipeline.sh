#!/bin/bash

echo ========================================================
echo ===			generating aneurysms 				  ===
echo ========================================================

vmtkcenterlines -ifile test_environment/SKD0050_baseline_model.vtp --pipe vmtkcenterlineattributes --pipe vmtkbranchextractor -ofile RCA_cl.vtp
#vmtkcenterlines -ifile test_environment/SKD0050_baseline_model.vtp -seedselector openprofiles -ofile total_cl.vtp
# vmtkcenterlines -ifile test_environment/SKD0050_baseline_model.vtp  -ofile total_cl.vtp


#vmtksurfacereader -ifile test_environment/SKD0050_baseline_model.vtp --pipe vmtkcenterlines --pipe vmtkrenderer --pipe vmtksurfaceviewer -opacity 0.25 --pipe vmtksurfaceviewer -i @vmtkcenterlines.o -array MaximumInscribedSphereRadius
