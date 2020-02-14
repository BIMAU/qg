#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
        mexErrMsgIdAndTxt("QG_mass:nrhs", "One input argument required.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("QG_mass:nlhs", "One output argument required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_mass:prhs", "Input must be a real uint64 scalar.");
    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));

    int ndim= (int)*mxGetPr(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)ndim, (mwSize)1, mxREAL);
    double *m = mxGetPr(plhs[0]);

    qg->mass(m);
}
