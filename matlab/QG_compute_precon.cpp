#include "mex.h"

#include <cstdint>

#include "QG.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgIdAndTxt("QG_compute_precon:nrhs", "One input argument required.");
    if (nlhs != 0)
        mexErrMsgIdAndTxt("QG_compute_precon:nlhs", "Zero output argument required.");

    // Get back the QG object
    if (mxGetNumberOfElements(prhs[0]) != 1 || mxGetClassID(prhs[0]) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("QG_compute_precon:prhs", "Input must be a real uint64 scalar.");

    QG::QG *qg = reinterpret_cast<QG::QG *>(*((uint64_t *)mxGetData(prhs[0])));

    qg->compute_precon();
}
