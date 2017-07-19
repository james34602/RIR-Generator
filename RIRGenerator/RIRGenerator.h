#include <math.h>
#include <stdlib.h>
#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#define M_PI 3.14159265358979323846
double Sinc(double x);
double microphoneType(double x, double y, double z, double* angle, int mtype);
void RIRCleanUp(double** imp, double** rr, int nMicrophones, double* LPI);
void rir_generator(double** imp, double c, int fs, int nMicrophones, double** rr, double* ss, double* LL, int betaNum, double* beta, int timeWidth, double* LPI, int nSamples, int microphone_type, int nOrder, int nDimension, double* angle, int isHighPassFilter);
double** gen2DArrayCALLOC(int arraySizeX, int arraySizeY);