# magnetometer_calib

Experimentation with magnetometer (hard iron) calibration. Collects samples of a 3D magnetometer 
and does a geometric sphere fit to find the center for bias/zero offset calibration of the axes.

Multiple samples are collected in a buffer with the condition that the new value must
keep a minimum distance to all previous ones. This forces the samples to occupy a large 
enough area on the samples sphere, so that the points are dispersed not only on a plane.

Sphere fitting is based on the algorithm shown in https://arxiv.org/ftp/arxiv/papers/1506/1506.02776.pdf. The
algorithm is modified so that in the first phase the sums are accumulated one by one, and 
finally the rest of the calculations are made as a final step. This way the samples don't need to be stored.

## Usage

 * Connect a magnetic sensor to an Arduino compatible processor and run the code.
 * Rotate the sensor around and samples will be printed to the terminal
 * After enough samples are collected, sphere fitting results are printed

You can also copy and paste the sample values to data.csv and run magnetometer_calib.m in Octave.
This will visualize the samples in a 3D plot.