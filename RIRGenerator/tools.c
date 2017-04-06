#include <math.h>
#include "sndfile.h"
void normalize(double* buffer, int num_samps, double maxval)
{
	double loudest_sample = 0.0;
	double multiplier = 0.0;
	int i;

	for (i = 0; i < num_samps; i++)
	{
		if (fabs(buffer[i]) > loudest_sample) loudest_sample = buffer[i];
	}

	multiplier = maxval / loudest_sample;

	for (i = 0; i < num_samps; i++)
	{
		buffer[i] *= multiplier;
	}
}
void channel_split(double* buffer, int num_frames, double** chan_buffers, int num_channels)
{
	int i;
	int samples = num_frames * num_channels;
	for (i = 0; i < samples; i++)
	{
		chan_buffers[(i % num_channels)][i / num_channels] = buffer[i];
	}
}

void channel_join(double** chan_buffers, int num_channels, double* buffer, int num_frames)
{
	int i;
	int samples = num_frames * num_channels;
	for (i = 0; i < samples; i++)
	{
		buffer[i] = chan_buffers[i % num_channels][i / num_channels];
	}
}
void make_mono(double* inbuffer, double* outbuffer, int num_channels, int num_frames)
{
	int i;
	for (i = 0; i < num_frames; i++)
	{
		int j;
		double sum = 0.0;
		for (j = 0; j < num_channels; j++)
		{
			sum += inbuffer[(i * num_channels) + j];
		}
		sum /= (double)num_channels;
		outbuffer[i] = sum;
	}
}
/* String tools */
int strcmpX(const char *str1, const char *str2)
{
	const unsigned char *ptr1 = (const unsigned char *)str1;
	const unsigned char *ptr2 = (const unsigned char *)str2;
	while (*ptr1 && *ptr1 == *ptr2) ++ptr1, ++ptr2;
	return (*ptr1 > *ptr2) - (*ptr2  > *ptr1);
}
char* strcatX(char *dest, const char *src)
{
	int len_1 = 0;
	int len_2 = 0;
	for (; dest[len_1] != '\0'; len_1++);
	for (; src[len_2] != '\0'; len_2++)
		dest[len_1 + len_2] = src[len_2];
	dest[len_1 + len_2] = '\0';
	return dest;
}