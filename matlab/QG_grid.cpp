#include "mex.h"
#include <cstdint>
#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgIdAndTxt("QG_grid:nrhs", "One input arguments required.");
    if (nlhs != 2)
        mexErrMsgIdAndTxt("QG_grid:nlhs", "Two output argument required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_grid:prhs", "Input must be a real uint64 scalar.");
    
    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));
    
    int m = qg->m() * qg->n();
    int n = 1;

    plhs[0] = mxCreateDoubleMatrix((mwSize)m, (mwSize)n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)m, (mwSize)n, mxREAL);
    
    double *x = mxGetPr(plhs[0]);
    double *y = mxGetPr(plhs[1]);

    qg->grid(x, y);
}
