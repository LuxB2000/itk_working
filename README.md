# ITK tools - DCE MRI pipelines
Some tools and bash pipelines to run registrations, segmentations and image manipulations with ITK. Tools have been designed in the DCE-MRI context. 

The tools are small binaries. Use the bash stript to run a pipeline. The main goal of these tools is to measure the uptake of contrast based agent in region of interests.

Warning: the framework has been designed for Linux stations. If you want to use it on MacOs or Window you need to write your own pipeline file. However, you should be able to compile them without error - not been tested so far.

## Requirements
CMake, ITK, VTK.

## Build
On Linux workstation:

> $ mkdir build bin

> $ cd build

> $ cmake ..

> $ make

## Tools
All the tools have a help, just run it without any input to have information. Most of the tools have been developped in the DCE-MRI contexts. It extracts region of interest, it does registration between volumes and volume and atlas. Some tools are dedicated to read and write images. Volumes are assumed to 3D but 4D volumes can be written and read. 4D volumes are timeline DCE-MRI, on each pixel location a vector contains values representing values along the time.

## Pipelines
