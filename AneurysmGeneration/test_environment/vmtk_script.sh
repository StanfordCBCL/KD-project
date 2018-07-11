#!/bin/bash 

echo 'hello!'

results_path=$1
base_path=${results_path%.vtp}
echo 'Using vmtk to analyze: ' $results_path

# vmtk will compute the centerline from scratch. 
# The default settings here employ an interactive point interactor that requires you to select points from the surface visualization that pops up as the seeds for the centerline. 
# This has been really finicky for me, and it's been hard for me to make sure that I'm actually selecting points on the surface (apparently, because I keep getting errors). 
# I'd probably recommend doing this a single centerline at a time to make sure you can just go through this once. 


######        if not computing from scratch, set the CL path here; otherwise, uncomment the lines below
cl_path=RCA_cl.vtp

#cl_path=$base_path'_cl.vtp'
#echo 'Output centerline will be stored at: ' $cl_path
#vmtkcenterlines -ifile $results_path --pipe vmtkcenterlineattributes --pipe vmtkbranchextractor -ofile $cl_path

# the next step is computing the reference system, which is basically just position along the centerline and the TNB coordinate system about the cetnerline path. 
ref_system_path=$base_path'_cl_rs.vtp'
echo 'Output reference system will be stored at: ' $ref_system_path
vmtkbifurcationreferencesystems -ifile $cl_path -radiusarray MaximumInscribedSphereRadius -blankingarray Blanking -groupidsarray GroupIds -ofile $ref_system_path

# this part splits the surface into any branches. 
clipped_path=$base_path'_clipped.vtp'
echo 'Branch clipping results will be stored at" ' $clipped_path
vmtkbranchclipper -ifile $results_path -centerlinesfile $cl_path -groupidsarray GroupIds -radiusarray MaximumInscribedSphereRadius -blankingarray Blanking -ofile $clipped_path

# This step actually performs the mapping of surface into the computed reference system. 
clipped_metrics=$base_path'_clipped_metrics.vtp'
echo 'Clipping metrics will be stored at: ' $clipped_metrics
vmtkbranchmetrics -ifile $clipped_path -centerlinesfile $cl_path -abscissasarray Abscissas -normalsarray ParallelTransportNormals -groupidsarray GroupIds -centerlineidsarray CenterlineIds -tractidsarray TractIds -blankingarray Blanking -radiusarray MaximumInscribedSphereRadius -ofile $clipped_metrics


# now, it's time to map the computed metrics from our results file into the right region with some corrections based on bifurcation regions. 
clipped_mapping=$base_path'_clipped_mapping.vtp'
echo 'Fine-tuned mapping results wil be applied to branches at: ' $clipped_mapping
vmtkbranchmapping -ifile $clipped_metrics -centerlinesfile $cl_path -referencesystemsfile $ref_system_path -normalsarray ParallelTransportNormals -abscissasarray Abscissas -groupidsarray GroupIds -centerlineidsarray CenterlineIds -tractidsarray TractIds -referencesystemsnormalarray Normal -radiusarray MaximumInscribedSphereRadius -blankingarray Blanking -angularmetricarray AngularMetric -abscissametricarray AbscissaMetric -ofile $clipped_mapping


# now is finally time to flatten the surface for visualization. 
# There are some parameters that we can tune: 
clipped_patching=$base_path'_clipped_patching.vtp'
patch_size=0.5
num_circ_patch=12
output_img='flattening.png'
echo 'Patching vtp is here: ' $clipped_patching
echo 'The longitudinal patchsize in mm: ' $patch_size
echo 'The number of subdivisions of -pi to pi: ' $num_circ_patch
echo 'The output image will be at: ' $output_img
vmtkbranchpatching -ifile $clipped_mapping -groupidsarray GroupIds -longitudinalmappingarray StretchedMapping -circularmappingarray AngularMetric -longitudinalpatchsize $patch_size -circularpatches $num_circ_patch -ofile $clipped_patching -patcheddatafile $output_img

