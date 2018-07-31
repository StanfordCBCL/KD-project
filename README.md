
## AneurysmGeneration

We want to artificially generate realistic aneurysms. 

## Current progress: 

The full pathway is 
1. generate aneurysms 
2. mesh the aneurysm model 
3. run simulation on Sherlock
4. post process results
5. map mesh and model points 
6. use mapping to extract relevant region and clip out the aneurysm for analysis 

`Aneurysm Generation` corresponds to step 1 , `batch_clip.sh` and `prelim_analysis.py` correspond to step 5-6. 


The desired analysis is: 
* Time-averaged Wall Shear Stress (TAWSS) as a function of Z-score/position
* Range of TAWSS as a function of Z-score/position
* CDF of TAWSS/area. 

## Methods 
### Dependencies: 
* `vtk`
* `numpy`
* `scipy`
* `matplotlib`
* `argparse`

### Generation 
A script has been written to manipulate the points of the model vtp file to produce an aneurysm-like shape which is radially symmetric about the cross section. It is possible to specify the aneurysm shape and length to fit different shape indices. It is possible to control the position along the vessel. 

The current method involves interpolation of radius as a function of path length, and does not consider path curvature or the underlying radial distribution in the Frenet-Serre frame; the result is potentially off-center aneurysms with non-smooth expansion and contraction at the inflow and outflow. Code has also been developed for 2-D interpolation that maps (s, theta) -> r, but the interpolating function is not well behaved. The current workaround is boolean expansion as max(original radius, new radius) coupled with smoothing in Simvascular after aneurysm generation. 

### Post Processing
Read in the parameters that were used to generate the actual aneurysm, then map the points from the mesh to the model to extract the right branch. From the right branch, we can isolate the points corresponding to the aneurysm. From this, we can define two `vtkPlane`s and cut the aneurysm. The true aneurysmal region, and not any other part of the model, will be identified using a `vtkConnectivityFilter`. 


### Todo: 
* The left sided centerline is questionable? And also the branching over there is much more complicated to handle. Will potentially develop a better method for introducing the aneurysm. 
* A current issue is that the smoothing operations take significant manual effort because the expansion does not produce new surface points, and thus stretches the existing elements. Potentially requires coupling with some sort of vtkdensifypoints or whatever it's called. 


