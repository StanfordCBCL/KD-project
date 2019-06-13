#!/bin/bash

# cp d1/int_r.dat ints_volumes 
# cp d1/Volume.dat ints_volumes
# mv ints_volumes/int_r.dat ints_volumes/int_r_d1.dat 
# mv ints_volumes/Volume.dat ints_volumes/Volume_d1.dat

# cp d2/int_r.dat ints_volumes 
# cp d2/Volume.dat ints_volumes
# mv ints_volumes/int_r.dat ints_volumes/int_r_d2.dat 
# mv ints_volumes/Volume.dat ints_volumes/Volume_d2.dat

# cp d3/int_r.dat ints_volumes 
# cp d3/Volume.dat ints_volumes
# mv ints_volumes/int_r.dat ints_volumes/int_r_d3.dat 
# mv ints_volumes/Volume.dat ints_volumes/Volume_d3.dat

# cp d4/int_r.dat ints_volumes 
# cp d4/Volume.dat ints_volumes
# mv ints_volumes/int_r.dat ints_volumes/int_r_d4.dat 
# mv ints_volumes/Volume.dat ints_volumes/Volume_d4.dat

# cp d5/int_r.dat ints_volumes 
# cp d5/Volume.dat ints_volumes
# mv ints_volumes/int_r.dat ints_volumes/int_r_d5.dat 
# mv ints_volumes/Volume.dat ints_volumes/Volume_d5.dat


# copy the baseline
scp Artificial/baseline/all_results.vtu Artificial/baseline/SKD0050_baseline_model-mesh-complete/mesh-complete.mesh.vtu Artificial/baseline/SKD0050_baseline_model-mesh-complete/mesh-surfaces/inflow.vtp alu2@login.sherlock.stanford.edu:/scratch/groups/amarsden/alex/baseline

# mesh-complete suffix
mc="mesh-complete.mesh.vtu"

# inflow suffix
inf="mesh-surfaces/inflow.vtp"

# target path
target="alu2@login.sherlock.stanford.edu:/scratch/groups/amarsden/alex"

# # left side 
# for shape in 2 4 6
# do 
# 	echo "using asi = $shape"
# 	for z in 1 2 3 4 5 
# 	do
# 		echo "using z = $z"

# 		results="Artificial/LAD/ASI$shape/lad$z/all_results.vtu"
# 		source="Artificial/LAD/ASI$shape/lad$z/SKD0050_baseline_model_modified_lad$z-mesh-complete/"

# 		echo $results
# 		echo $source
# 		scp $results "$source$mc" "$source$inf" "$target/lad/asi$shape/lad$z"
# 	done
# done 


right side 
for shape in 4 6
do 
	echo "using asi = $shape"
	for z in 1 2 3 4 5 
	do
		echo "using z = $z"

		results="Artificial/RCA/ASI$shape/p$z/all_results.vtu"
		source="Artificial/RCA/ASI$shape/p$z/SKD0050_baseline_model_modified_p$z-mesh-complete/"

		echo $results
		echo $source
		scp $results "$source$mc" "$source$inf" "$target/rca/asi$shape/p$z"
	done
done 

# for shape in 2
# do 
# 	echo "using asi = $shape"
# 	for z in 1 2 3 4 5 
# 	do
# 		echo "using z = $z"

# 		for pos in p m d 
# 		do
# 			results="Artificial/RCA/ASI$shape/$pos$z/all_results.vtu"
# 			source="Artificial/RCA/ASI$shape/$pos$z/SKD0050_baseline_model_modified_$pos$z-mesh-complete/"

# 			echo $results
# 			echo $source
# 			scp $results "$source$mc" "$source$inf" "$target/rca/asi$shape/$pos$z"
# 		done
# 	done
# done 


