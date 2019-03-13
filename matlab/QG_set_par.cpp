#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3)
        mexErrMsgIdAndTxt("QG_set_par:nrhs", "Three input arguments required.");
    if (nlhs != 0)
        mexErrMsgIdAndTxt("QG_set_par:nlhs", "Zero output arguments required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_set_par:prhs", "Input must be a real uint64 scalar.");
    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));

    int par = (int)*mxGetPr(prhs[1]);
    double val = (double)*mxGetPr(prhs[2]);

    qg->set_par(par, val);
}
