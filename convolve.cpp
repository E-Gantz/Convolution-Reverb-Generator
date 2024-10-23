#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <string>

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE   16

/*  Standard sample size in bytes  */		
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC        1
#define STEREOPHONIC      2

#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

//NOTE: according to http://soundfile.sapp.org/doc/WaveFormat/ " 16-bit samples are stored as 2's-complement signed integers, ranging from -32768 to 32767. "

//convert a 16-bit 2's-complement signed integer to a double scaled to -1.0 to +1.0
double sampleTodouble(int16_t value){
    return double(value) / (double)32767.0;
}

//convert a double scaled to -1.0 to +1.0 to a 16-bit 2's-complement signed integer
int16_t doubleToSample(double value){
    int16_t maxVal = 32767;
    int16_t result;
    if (value >= 1.0){
        result = maxVal;
    }
    else if (value <= -1.0){
        result = -maxVal;
    }
    else {
        result = rint((double)value * (double)maxVal);
    }

    return result;
}

//convert the vector of wav samples to doubles
std::vector<double> samplesTodoubles(std::vector<int16_t> samples){
   std::vector<double> doubleSamples;
   for (int16_t sample : samples){
        double convertedSample = sampleTodouble(sample);
        doubleSamples.push_back(convertedSample);
   }
   return doubleSamples;
}

//convert the vector of doubles into wav samples
std::vector<int16_t> doublesToSamples(std::vector<double> doubleSamples){
   std::vector<int16_t> samples;
   for (double sample : doubleSamples){
        int16_t convertedSample = doubleToSample(sample);
        samples.push_back(convertedSample);
    }
   return samples;
}


std::vector<int16_t> readWavFile(char *filename){
    std::ifstream file(filename, std::ios::binary);

    if (!file.is_open()){
        fprintf(stdout, "File %s cannot be opened for reading\n", filename);
        return std::vector<int16_t>();
    }

    //read header
    char header[44];
    file.read(header, 44);

    //Check if the file is mono, 16-bit, 44.1 kHz
    //header[22] is numchannels, 34 is bits per sample, 24 and 26 is sample rate
    if (header[22] != 1 || header[34] != 16 || *reinterpret_cast<int*>(header + 24) != 44100){
        std::cerr << "Unsupported WAV file format. Must be mono, 16-bit, 44.1 kHz." << std::endl;
        return std::vector<int16_t>();
    }

    //store samples
    std::vector<int16_t> samples;
    int16_t sample;
    while (file.read(reinterpret_cast<char*>(&sample), sizeof(int16_t))){
        samples.push_back(sample);
    }

    file.close();
    return samples;
}


size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}

size_t fwriteShortLSB(int16_t data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}

void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    int16_t frameSize = channels * BYTES_PER_SAMPLE;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((int16_t)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}

//write audio data to a new WAV file
void writeWavFile(char *filename, std::vector<int16_t> samples, int numberOfChannels){
    FILE *outputFileStream = fopen(filename, "wb");
    if (outputFileStream == NULL){
        fprintf(stdout, "File %s cannot be opened for writing\n", filename);
        return;
    }

    int numberOfSamples = samples.size();

    /*  Write the WAVE file header  */
    writeWaveFileHeader(numberOfChannels, numberOfSamples,
                        SAMPLE_RATE, outputFileStream);

    //Write audio data (samples) to the file
    for (int i = 0; i < numberOfSamples; i++){
        fwriteShortLSB(samples[i], outputFileStream);
    }

    fclose(outputFileStream);
}

//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(double data[], unsigned long nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

int nearestPower(int M){
    int i = 0;
    while (true){
        if (M <= pow(2, i)){
            return pow(2,i);
        }
        else {
            i++;
        }
    }
}


void convolve(std::vector<double> input, std::vector<double> filter, std::vector<double> &y){
    int N = input.size();
    int M = filter.size();
    int m = M-1;
    int output_size = (M);
    int paddingSize = nearestPower(output_size);
    int bigPad = 2*paddingSize;
    int numChunkSamples = (paddingSize-(m));
    int numChunks = N/numChunkSamples;
    std::vector<double> OLAP(m, 0.0);

    double* paddedInput = new double[bigPad];
    double* paddedImpulse = new double[bigPad];
    double* paddedOutput = new double[paddingSize];

    for(int i=0; i<bigPad; i++){
        paddedImpulse[i] = 0.0;
    }
    int count = 0;
    for (int i=0; i<M; i++){
        paddedImpulse[count] = filter[i];
        count += 2;
    }
    four1(paddedImpulse-1, paddingSize, 1);
    int inputPos = 0;
    int outputPos = 0;
    for (int i=0; i<numChunks; i++){
        for(int j=0; j<bigPad; j++){
            paddedInput[j] = 0.0;
        }
        int counter = 0;
        for (int j=0; j<paddingSize-(m); j++){
            paddedInput[counter] = input[inputPos];
            inputPos++;
            counter += 2;
        }
        four1(paddedInput-1, paddingSize, 1);
        
        //complex multiplication
        for (int j=0; j<bigPad; j+=2){
            double temp = (paddedInput[j]*paddedImpulse[j]) - (paddedInput[j+1]*paddedImpulse[j+1]);
            paddedInput[j+1] = (paddedInput[j]*paddedImpulse[j+1]) + (paddedInput[j+1]*paddedImpulse[j]);
            paddedInput[j] = temp;
        }
        four1(paddedInput-1, paddingSize, -1);
        counter = 0;
        for (int j=0; j<paddingSize; j++){
            paddedOutput[j] = paddedInput[counter];
            paddedOutput[j] = (double)paddedOutput[j] / (double)N;
            paddedOutput[j] = (double)paddedOutput[j] * (double)1.7;    //turn the gain up a little bit
            counter +=2;
        }
        for (int j=0; j<m; j++){
            paddedOutput[j] = paddedOutput[j] + OLAP[j];
        }
        for (int j = numChunkSamples; j<numChunkSamples+(m); j++){
            OLAP[j-numChunkSamples] = paddedOutput[j];
        }
        for (int j=0; j<numChunkSamples; j++){
            y[outputPos] = paddedOutput[j];
            outputPos++;
        }
    }
    for (int j=0; j<m; j++){
        y[outputPos] = OLAP[j];
        outputPos++;
    }
}

int main(int argc, char *argv[]){
    char *inputFile = NULL;
    char *IRFile = NULL;
    char *outputFile = NULL;

    if (argc != 4){
        std::cout << "Usage: " << argv[0] << " inputfile IRfile outputfile" << std::endl;
        return 1;
    }

    inputFile = argv[1];
    IRFile = argv[2];
    outputFile = argv[3];

    //get the data section of the input file as a vector
    std::vector<int16_t> wavSamples = readWavFile(inputFile);
    if (wavSamples.empty()){
        return 1;
    }
    std::vector<int16_t> irSamples = readWavFile(IRFile);
    if (irSamples.empty()){
        return 1;
    }

    //ints to doubles
    std::vector<double> doubleSamples = samplesTodoubles(wavSamples);
    std::vector<double> doubleIR = samplesTodoubles(irSamples);
    std::vector<double> convolved(doubleSamples.size() + doubleIR.size() - 1, 0.0);

    //convolution
    convolve(doubleSamples, doubleIR, convolved);

    //doubles back to ints
    std::vector<int16_t> outputSamples = doublesToSamples(convolved);

    //Write the converted samples to a new WAV file
    writeWavFile(outputFile, outputSamples, 1);

    return 0;
}