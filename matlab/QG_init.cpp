#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
        mexErrMsgIdAndTxt("QG_init:nrhs", "Two input arguments required.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("QG_init:nlhs", "One output argument required.");

    int nx = (int)*mxGetPr(prhs[0]);
    int ny = (int)*mxGetPr(prhs[1]);

    QG::QG *qg = new QG::QG(nx, ny);

    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *((uint64_t *)mxGetData(plhs[0])) = reinterpret_cast<uint64_t>(qg);
}
