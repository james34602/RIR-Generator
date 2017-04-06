#include <stdio.h>
#include "RIRGenerator.h"
#include "tools.h"
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
		"               and MicrophonesXXImpulseResponse.wav\n");
	SNDFILE* wavOut = NULL;
	SF_INFO impProp;
	int i, j, k;
	double c;
	int SampleRate;
	int nMicrophones;
	double** imp;
	double** rr;
	double* beta;
	double *sndBuf;
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
	int totalFrames;
	double fillrr;
	printf("Sound velocity:");
	scanf("%lf", &c);
	printf("Sample rate:");
	scanf("%d", &SampleRate);
	printf("Numbers of microphones:");
	scanf("%d", &nMicrophones);
	rr = gen2DArrayCALLOC(nMicrophones, 3);
	for (i = 0; i < nMicrophones; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if(j==0)
				printf("Fill in x position of microphone %d: ", i + 1);
			else if (j == 1)
				printf("Fill in y position of microphone %d: ", i + 1);
			else if (j == 2)
				printf("Fill in z position of microphone %d: ", i + 1);
			scanf("%lf", &fillrr);
			rr[i][j] = fillrr;
		}
	}
	printf("Source x,y,z:");
	for (i = 0; i < 3; i++)
	{
		scanf("%lf", &fillrr);
		ss[i] = fillrr;
	}
	printf("Room size x,y,z:");
	for (i = 0; i < 3; i++)
	{
		scanf("%lf", &fillrr);
		LL[i] = fillrr;
	}
	printf("Numbers of beta:");
	scanf("%d", &betaNum);
	beta = gen1DArrayCALLOC(betaNum);
	for (i = 0; i < betaNum; i++)
	{
		printf("Beta:");
		scanf("%lf", &fillrr);
		beta[i] = fillrr;
	}
	printf("Total samples you want to generate:");
	scanf("%d", &nSamples);
	printf("Microphone types[Bidirectional=0, hypercardioid=1, cardioid=2, subcardioid=3, omnidirectional=4]:");
	scanf("%d", &microphone_type);
	printf("Reflection order, -1 is maximum order:");
	scanf("%d", &nOrder);
	printf("Room dimensions, 2=2D 3=3D:");
	scanf("%d", &nDimension);
	if (nDimension < 2)
		nDimension = 2;
	else
		nDimension = 3;
	printf("Enable orientations? 0=Disable 1=Enable:");
	scanf("%d", &orientationEnable);
	if (orientationEnable == 1)
	{
		printf("Orientation:");
		for (i = 0; i < 2; i++)
		{
			scanf("%lf", &fillrr);
			orientation[i] = fillrr;
		}
	}
	printf("Need high pass filtering of the impulse response? 0=Disable 1=Enable:");
	scanf("%d", &isHighPassFilter);
	imp = rir_generator(c, SampleRate, nMicrophones, rr, ss, LL, betaNum, beta, nSamples, microphone_type, nOrder, nDimension, orientationEnable, orientation, isHighPassFilter);
	char dynamicName[256];
	sprintf(dynamicName, "Microphone%iImpulseResponse.csv", nMicrophones);
	FILE *impFile = fopen(dynamicName, "w+");
	for (i = 0; i < nSamples; i++)
	{
		for (k = 0; k < nMicrophones; k++)
		{
			if (k == nMicrophones - 1) {
				if(imp[k][i] == 0.0)
					fprintf(impFile, "0");
				else
					fprintf(impFile, "%4.65f", imp[k][i]);
			}
			else {
				if (imp[k][i] == 0.0)
					fprintf(impFile, "0,");
				else
					fprintf(impFile, "%4.65f,", imp[k][i]);
			}
		}
		fprintf(impFile, "\n");
	}
	fclose(impFile);
	sprintf(dynamicName, "Microphones%iImpulseResponse.wav", nMicrophones);
	totalFrames = nMicrophones * nSamples;
	sndBuf = (double*)calloc(totalFrames, sizeof(double));
	channel_join(imp, nMicrophones, sndBuf, nSamples);
	normalize(sndBuf, totalFrames, 0.9);
	impProp.samplerate = SampleRate;
	impProp.channels = nMicrophones;
	impProp.frames = totalFrames;
	impProp.format = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
	wavOut = sf_open(dynamicName, SFM_WRITE, &impProp);
	sf_writef_double(wavOut, sndBuf, nSamples);
	sf_close(wavOut);
	free(beta);
	for (int i = 0; i < nMicrophones; i++)
		free(rr[i]);
	free(rr);
	for (int i = 0; i < nMicrophones; i++)
		free(imp[i]);
	free(imp);
	free(sndBuf);
	exit(0);
}