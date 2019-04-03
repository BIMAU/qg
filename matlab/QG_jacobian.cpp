#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3)
        mexErrMsgIdAndTxt("QG_jacobian:nrhs", "Two input arguments required.");
    if (nlhs != 3)
        mexErrMsgIdAndTxt("QG_jacobian:nlhs", "Three output arguments required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_jacobian:prhs", "Input must be a real uint64 scalar.");
    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));

    double *x = mxGetPr(prhs[1]);
    int mp1 = mxGetM(prhs[1]) + 1;
    double sig= (double)*mxGetPr(prhs[2]);
   
    plhs[0] = mxCreateNumericMatrix(mp1, 1, mxINT32_CLASS, mxREAL);
    int *beg = (int *)mxGetPr(plhs[0]);
    int nco=4*(mp1-1)*9;
    plhs[1] = mxCreateNumericMatrix(nco, 1, mxINT32_CLASS, mxREAL);
    int *jco = (int *)mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix((mwSize)nco, (mwSize)1, mxREAL);
    double *co = mxGetPr(plhs[2]);

    qg->jacobian(x, sig, beg, jco, co);
}
