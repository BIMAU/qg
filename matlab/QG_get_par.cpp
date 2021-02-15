#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
        mexErrMsgIdAndTxt("QG_get_par:nrhs", "Two input arguments required.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("QG_get_par:nlhs", "One output argument required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_get_par:prhs", "Input must be a real uint64 scalar.");
    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));

    int par = (int)*mxGetPr(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    double *val = mxGetPr(plhs[0]);

    *val = qg->get_par(par);
}
