#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3)
        mexErrMsgIdAndTxt("QG_init:nrhs", "Wrong number of input arguments.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("QG_init:nlhs", "One output argument required.");

    int nx = (int)*mxGetPr(prhs[0]);
    int ny = (int)*mxGetPr(prhs[1]);
    int pr = (int)*mxGetPr(prhs[2]);

    QG::QG *qg = new QG::QG(nx, ny, pr);

    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *((uint64_t *)mxGetData(plhs[0])) = reinterpret_cast<uint64_t>(qg);
}
