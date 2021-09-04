
#ifndef __INCLUDE_MYFIO__
#define __INCLUDE_MYFIO__
#include "myUtils.h"
#include <hdf5.h>


INT get_REAL8TimeSeries_waveform(REAL8TimeSeries ** hplus, REAL8TimeSeries ** hcross, FILE * fp);
INT get_REAL8TimeSeries(REAL8TimeSeries ** Series, FILE * fp);
INT get_COMPLEX16TimeSeries_waveform(COMPLEX16TimeSeries ** hout, FILE * fp);

INT cmd_mkdir(CHAR *folderName);
void CreateFolder(CHAR *folderName);
INT read_waveform(REAL8Vector **time, 
                  REAL8Vector **hreal, 
                  REAL8Vector **himag,
                  FILE *file);

INT DumpREAL8ArrayTohdf5(CHAR *fname, CHAR *dname, REAL8Array *array, INT is_delete);


#endif
