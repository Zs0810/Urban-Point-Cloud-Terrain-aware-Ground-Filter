# Urban-Point-Cloud-Terrain-aware-Ground-Filter
MATLAB implementation of terrain-aware ground filtering for urban photogrammetric point clouds.

## Description
This code implements an adaptive ground filtering method for urban photogrammetric point clouds.
The pipeline includes:
- Adaptive bottom-layer identification
- Progressive multi-scale and multi-direction morphological filtering
- Terrain-prior constrained ground point refinement
- DEM reconstruction from filtered ground points

## Requirements
- MATLAB (R2020b or newer)
- Image Processing Toolbox

## Usage
1. Prepare your point cloud file (`.txt` format: X, Y, Z)
2. Modify the input path parameters in the main script
3. Run the main script
4. Output: classified ground/non-ground points and generated DEM

## Input
Text file with 3 columns: X, Y, Z

## Output
- Ground point cloud
- Non-ground point cloud
- DEM raster

## Citation
If you use this code for your research, please cite our paper:
