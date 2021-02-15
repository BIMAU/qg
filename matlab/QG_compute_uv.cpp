#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
        mexErrMsgIdAndTxt("QG_compute_uv:nrhs", "Two input arguments required.");
    if (nlhs != 2)
        mexErrMsgIdAndTxt("QG_compute_uv:nlhs", "Two output arguments required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_compute_uv:prhs", "Input must be a real uint64 scalar.");
    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));

    double *x = mxGetPr(prhs[1]);
    int m     = mxGetM(prhs[1]);
    int n     = mxGetN(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)m/2, (mwSize)n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)m/2, (mwSize)n, mxREAL);
    double *u = mxGetPr(plhs[0]);
    double *v = mxGetPr(plhs[1]);

    qg->compute_uv(x, u, v);
}
