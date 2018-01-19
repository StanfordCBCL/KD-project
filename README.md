
## AneurysmGeneration

We want to artificially generate realistic aneurysms. 

## Current progress: 

A script has been written to expand the shell points to produce an aneurysm-like shape which is radially symmetric about the cross section. 

Further steps will be taken to generate different kinds of shapes. 

![simple expansion](AneurysmGeneration/screenshots/progress_1_12.png)

Progress was also made to increase the smoothness of the aneurysm surface. Initially, the number of centerline points was insufficient, and many points along the vessel surface were mapping to the same centerline coordinate. Increasing the number of centerline points allows for a smoother-looking surface post-expansion. 

![smoother expansion](AneurysmGeneration/screenshots/progress_1_12_2.png)

### Notes:
* discovered that the model looked ribbed because there were not enough points in the centerline; centerline points increased. 

### Next steps: 
* work on a real model 
* develop centerline processing to read in .path file
* update interpolation to force a smoother blend 
* develop processing for automated read of coronary tree and more robust positioning of aneurysm growth. 