#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3)
        mexErrMsgIdAndTxt("QG_bilin:nrhs", "Two input arguments required.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("QG_bilin:nlhs", "One output argument required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_bilin:prhs", "Input must be a real uint64 scalar.");
    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));


    double *x = mxGetPr(prhs[1]);
    int mx = mxGetM(prhs[1]);
    int nx = mxGetN(prhs[1]);
    
    double *y = mxGetPr(prhs[2]);
    //  int my = mxGetM(prhs[2]);
    // int ny = mxGetN(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)mx, (mwSize)nx, mxREAL);
    double *z = mxGetPr(plhs[0]);

    qg->bilin(x, y, z);
}
