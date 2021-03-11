#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgIdAndTxt("QG_get_par:nrhs", "One input argument required.");
    if (nlhs != 3)
        mexErrMsgIdAndTxt("QG_get_par:nlhs", "Three output arguments required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_get_nondim:prhs", "Input must be a real uint64 scalar.");
    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));


    plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    double *Lxdim = mxGetPr(plhs[0]);
    double *Lydim = mxGetPr(plhs[1]);
    double *Udim  = mxGetPr(plhs[2]);

    *Lxdim = qg->Lxdim();
    *Lydim = qg->Lydim();
    *Udim  = qg->Udim();
}
