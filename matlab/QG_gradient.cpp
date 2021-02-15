#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
        mexErrMsgIdAndTxt("QG_gradient:nrhs", "Two input arguments required.");
    if (nlhs != 2)
        mexErrMsgIdAndTxt("QG_gradient:nlhs", "Two output argument required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_gradient:prhs", "Input must be a real uint64 scalar.");
    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));

    double *field = mxGetPr(prhs[1]);
    int m = mxGetM(prhs[1]);
    int n = mxGetN(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)m, (mwSize)n, mxREAL);
    double *gradx = mxGetPr(plhs[0]);
    double *grady = mxGetPr(plhs[1]);

    qg->gradient(field, gradx, grady);
}
