#pragma once
#include "sndfile.h"
void normalize(double* buffer, int num_samps, double maxval);
void channel_split(double* buffer, int num_frames, double** chan_buffers, int num_channels);
void channel_join(double** chan_buffers, int num_channels, double* buffer, int num_frames);
void make_mono(double* inbuffer, double* outbuffer, int num_channels, int num_frames);
int strcmpX(const char *str1, const char *str2);
char* strcatX(char *dest, const char *src);