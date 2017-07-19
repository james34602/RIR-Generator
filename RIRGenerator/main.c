#include <stdio.h>
#include "RIRGenerator.h"
#define HAVESNDFILE
#ifdef HAVESNDFILE
#include "sndfile.h"
#endif
void normalise(double* buffer, int num_samps, double maxval)
{
	double loudest_sample = 0.0, multiplier = 0.0;
	int i;
	for (i = 0; i < num_samps; i++)
		if (fabs(buffer[i]) > loudest_sample) loudest_sample = buffer[i];
	multiplier = maxval / loudest_sample;
	for (i = 0; i < num_samps; i++)
		buffer[i] *= multiplier;
}
void channel_join(double** chan_buffers, int num_channels, double* buffer, int num_frames)
{
	int i, samples = num_frames * num_channels;
	for (i = 0; i < samples; i++)
		buffer[i] = chan_buffers[i % num_channels][i / num_channels];
}
int main()
{
	printf("--------------------------------------------------------------------\n"
		"| Room Impulse Response Generator                                  |\n"
		"|                                                                  |\n"
		"| Computes the response of an acoustic source to one or more       |\n"
		"| microphones in a reverberant room using the image method [1,2].  |\n"
		"|                                                                  |\n"
		"| Author    : dr.ir. Emanuel Habets (ehabets@dereverberation.org)  |\n"
		"| C implementation : James Fung(James34602) (james34602@gmail.com) |\n"
		"--------------------------------------------------------------------\n\n"
		"Input parameters:\n"
		" c           : sound velocity in m/s.\n"
		" SampleRate  : sampling frequency in Hz.\n"
		" r           : N x 3 array specifying the (x,y,z) coordinates of the\n"
		"               receiver(s) in m.\n"
		"               [x1 y1 z1 x2 y2 z2 ... xn yn zn]\n"
		" s           : 1 x 3 vector specifying the (x,y,z) coordinates of the\n"
		"               source in m.\n"
		" L           : 1 x 3 vector specifying the room dimensions (x,y,z) in m.\n"
		" beta        : 1 x 6 vector specifying the reflection coefficients\n"
		"               [beta_x1 beta_x2 beta_y1 beta_y2 beta_z1 beta_z2 ... beta_zn] or\n"
		"               beta = reverberation time (T_60) in seconds.\n"
		" nsample     : number of samples to calculate.\n"
		" mtype       : omnidirectional, subcardioid, cardioid, hypercardioid,\n"
		"               bidirectional.\n"
		" order       : reflection order, i.e. maximum order = -1.\n"
		" dim         : room dimension (2 or 3)\n"
		" orientation : direction in which the microphones are pointed, specified using\n"
		"               azimuth and elevation angles (in radians), default is [0 0].\n"
		" hp_filter   : use 'false' to disable high-pass filter\n"
		"Output parameters:\n"
		" h           : M x nsample matrix containing the calculated room impulse\n"
		"               response(s).Output as MicrophoneXXImpulseResponse.csv\n"
#ifdef HAVESNDFILE
		"               and MicrophonesXXImpulseResponse.wav\n");
	SNDFILE* wavOut = NULL;
	SF_INFO impProp;
	double *sndBuf;
	int totalFrames;
#else
		);
#endif
	int i, j, k;
	double c;
	int SampleRate;
	int nMicrophones;
	double ss[3];
	double LL[3];
	int nSamples;
	int microphone_type;
	int nOrder;
	int nDimension;
	int betaNum;
	int orientationEnable;
	double orientation[2];
	int isHighPassFilter;
	double reverberation_time = 0;
	// Below should not be changed after initalise
	printf("Sound velocity:");
	scanf("%lf", &c);
	printf("Sample rate:");
	scanf("%d", &SampleRate);
	printf("Numbers of microphones:");
	scanf("%d", &nMicrophones);
	printf("Room size x,y,z:");
	for (i = 0; i < 3; i++)
		scanf("%lf", &LL[i]);
	double beta[6] = { 0,0,0,0,0,0 };
	printf("Numbers of beta:");
	scanf("%d", &betaNum);
	for (i = 0; i < betaNum; i++)
	{
		printf("Beta:");
		scanf("%lf", &beta[i]);
	}
	printf("Total samples you want to generate(if numbers of beta > 1, this can be set to zero for auto generation):");
	scanf("%d", &nSamples);
	if (!nSamples && betaNum > 1)
	{
		double V = LL[0] * LL[1] * LL[2];
		double alpha = ((1 - pow(beta[0], 2)) + (1 - pow(beta[1], 2)))*LL[1] * LL[2] +
			((1 - pow(beta[2], 2)) + (1 - pow(beta[3], 2)))*LL[0] * LL[2] +
			((1 - pow(beta[4], 2)) + (1 - pow(beta[5], 2)))*LL[0] * LL[1];
		reverberation_time = 24 * log(10.0)*V / (c*alpha);
		if (reverberation_time < 0.128)
			reverberation_time = 0.128;
		nSamples = (int)(reverberation_time * SampleRate);
	}
	else if (!nSamples && betaNum <= 1)
	{
		printf("Cannot generate zero sample impulse response, at least numbers of beta should larger than 1");
		return 0;
	}
	if (betaNum == 1)
	{
		double V = LL[0] * LL[1] * LL[2];
		double S = 2 * (LL[0] * LL[2] + LL[1] * LL[2] + LL[0] * LL[1]);
		reverberation_time = beta[0];
		if (reverberation_time != 0)
		{
			double alfa = 24 * V*log(10.0) / (c*S*reverberation_time);
			if (alfa > 1)
				printf("Error: The reflection coefficients cannot be calculated using the current room parameters, i.e. room size and reverberation time.\nPlease "
					"specify the reflection coefficients or change the room parameters.");
			for (int i = 0; i < 6; i++)
				beta[i] = sqrt(1 - alfa);
		}
		else
		{
			for (int i = 0; i < 6; i++)
				beta[i] = 0;
		}
	}
	double** rr = gen2DArrayCALLOC(nMicrophones, 3);
	int timeWidth = 2 * ROUND(0.004*SampleRate); // The width of the low-pass FIR equals 8 ms;
	double* LPI = (double*)calloc(timeWidth, sizeof(double));
	double** imp = gen2DArrayCALLOC(nMicrophones, nSamples);
	// Below parameters can change whatever you want
	printf("Sound source x,y,z:");
	for (i = 0; i < 3; i++)
		scanf("%lf", &ss[i]);
	for (i = 0; i < nMicrophones; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (j == 0)
				printf("Fill in x position of microphone %d: ", i + 1);
			else if (j == 1)
				printf("Fill in y position of microphone %d: ", i + 1);
			else if (j == 2)
				printf("Fill in z position of microphone %d: ", i + 1);
			scanf("%lf", &rr[i][j]);
		}
	}
	printf("Microphone types[Bidirectional=0, hypercardioid=1, cardioid=2, subcardioid=3, omnidirectional=4]:");
	scanf("%d", &microphone_type);
	printf("Room dimensions, 2=2D 3=3D:");
	scanf("%d", &nDimension);
	if (nDimension == 2)
	{
		beta[4] = 0;
		beta[5] = 0;
	}
	else
		nDimension = 3;
	printf("Enable orientations? 0=Disable 1=Enable:");
	scanf("%d", &orientationEnable);
	if (orientationEnable == 1)
	{
		printf("Orientation 1:");
		scanf("%lf", &orientation[0]);
		printf("Orientation 2:");
		scanf("%lf", &orientation[1]);
	}
	else
	{
		orientation[0] = 0;
		orientation[1] = 0;
	}
	printf("Reflection order, -1 is maximum order:");
	scanf("%d", &nOrder);
	if (nOrder < -1)
		nOrder = -1;
	printf("High pass filtering of the impulse response? 0=Disable 1=Enable:");
	scanf("%d", &isHighPassFilter);
	rir_generator(imp, c, SampleRate, nMicrophones, rr, ss, LL, betaNum, beta, timeWidth, LPI, nSamples, microphone_type, nOrder, nDimension, orientation, isHighPassFilter);
	char dynamicName[256];
	sprintf(dynamicName, "Microphone%iImpulseResponse.csv", nMicrophones);
	FILE *impFile = fopen(dynamicName, "w+");
	for (i = 0; i < nSamples; i++)
	{
		for (k = 0; k < nMicrophones; k++)
		{
			if (k == nMicrophones - 1) {
				if (imp[k][i] == 0.0)
					fprintf(impFile, "0");
				else
					fprintf(impFile, "%4.18f", imp[k][i]);
			}
			else {
				if (imp[k][i] == 0.0)
					fprintf(impFile, "0,");
				else
					fprintf(impFile, "%4.18f,", imp[k][i]);
			}
		}
		fprintf(impFile, "\n");
	}
	fclose(impFile);
#ifdef HAVESNDFILE
	sprintf(dynamicName, "Microphones%iImpulseResponse.wav", nMicrophones);
	totalFrames = nMicrophones * nSamples;
	sndBuf = (double*)calloc(totalFrames, sizeof(double));
	channel_join(imp, nMicrophones, sndBuf, nSamples);
	normalise(sndBuf, totalFrames, 0.9);
	impProp.samplerate = SampleRate;
	impProp.channels = nMicrophones;
	impProp.frames = nSamples;
	impProp.format = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
	wavOut = sf_open(dynamicName, SFM_WRITE, &impProp);
	sf_writef_double(wavOut, sndBuf, nSamples);
	sf_close(wavOut);
	free(sndBuf);
#endif
	RIRCleanUp(imp, rr, nMicrophones, LPI);
	return 0;
}