#include <math.h>
#include <stdlib.h>
#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#define M_PI 3.14159265358979323846
double Sinc(double x);
double microphoneType(double x, double y, double z, double* angle, int mtype);
double** rir_generator(double c, int fs, int nMicrophones, double** rr, double* ss, double* LL, int betaNum, double* beta, int nSamples, int microphone_type, int nOrder, int nDimension, int orientationEnable, double* orientation, int isHighPassFilter);
double** gen2DArrayCALLOC(int arraySizeX, int arraySizeY);
double* gen1DArrayCALLOC(int arraySize);