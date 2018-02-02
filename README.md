
## AneurysmGeneration

We want to artificially generate realistic aneurysms. 

## Current progress: 

A script has been written to expand the shell points to produce an aneurysm-like shape which is radially symmetric about the cross section. 

Further steps will be taken to generate different kinds of shapes. 

![simple expansion](AneurysmGeneration/screenshots/progress_1_12.png)

Progress was also made to increase the smoothness of the aneurysm surface. Initially, the number of centerline points was insufficient, and many points along the vessel surface were mapping to the same centerline coordinate. Increasing the number of centerline points allows for a smoother-looking surface post-expansion. 

![smoother expansion](AneurysmGeneration/screenshots/progress_1_12_2.png)

Implementation of clamped boundary conditions allows for a smoother transition zone (zero first derivative) into healthy coronary vessel. 

![clamped bc example](AneurysmGeneration/screenshots/progress_1_19.png)

As of 2/1/2018: Full model generation kind of works! Shapes are awkward but full model is workable. 

![expansion of RCA](AneurysmGeneration/screenshots/progress_2_1.png)
![WT of RCA](AneurysmGeneration/screenshots/progress_2_1_2.png)




### Notes:
* discovered that the model looked ribbed because there were not enough points in the centerline; centerline points increased.
* ridge when bc is not 0 looks really awkward; now using 1D cubic spline interpolation allowing zero-derivative boundary conditions.  



### Next steps: 
* 