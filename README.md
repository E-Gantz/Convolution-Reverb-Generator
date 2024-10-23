# Convolution Reverb Generator
This is a convolution reverb program, it takes two audio files as input and produces a third audio file as output.

The first file should be whatever audio you want to add reverb to, the second file should be an Impulse Response file for whatever reverb you want to add.

The program then converts the samples from the input files (which are stored as integers) to doubles scaled to -1.0 to 1.0 and uses a Fast Fourier Transform, complex multiplication, and an inverse FFT to perform convolution on the two sets of samples to produce a third set, which are then converted back to integers and written to the output file.

Currently this is a command line program that only supports 16bit 44100Hz mono wav files

## Commands
To compile: `g++ -O3 -o convolve convolve.cpp`

To run: `./convolve <inputfile>.wav <IRFile>.wav <OutputFile>.wav`