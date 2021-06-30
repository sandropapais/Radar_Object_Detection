# Radar_Object_Detection
 
## FMCW Radar Waveform Design

Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp. For given system requirements the calculated slope should be around 2e13

## Target Generatio and Simulation

Simulate Target movement and calculate the beat or mixed signal for every timestep. A beat signal should be generated such that once range FFT implemented, it gives the correct range i.e the initial position of target assigned with an error margin of +/- 10 meters.

## Range FFT (1st FFT)

Implement the Range FFT on the Beat or Mixed Signal and plot the result. A correct implementation should generate a peak at the correct range, i.e the
initial position of target assigned with an error margin of +/- 10 meters.

## 2D CFAR

Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map. The 2D CFAR processing should be able to suppress the noise and separate
the target signal. The output should match the image shared in walkthrough.

Implementation steps for the 2D CFAR process.

Selection of Training, Guard cells and offset.

Steps taken to suppress the non-thresholded cells at the edges.
