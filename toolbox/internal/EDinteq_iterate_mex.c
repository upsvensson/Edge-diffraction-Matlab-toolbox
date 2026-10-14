#include "mex.h"
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

/*
 * EDinteq_iterate_mex.c
 *
 * Qfinal = EDinteq_iterate_mex(in_idx, out_idx, coeff, Q_firstterm, niterations, Nout)
 *
 * Inputs:
 *   in_idx       uint32 vector, MATLAB 1-based indices
 *   out_idx      uint32 vector, MATLAB 1-based indices
 *   coeff        double vector, real or complex
 *   Q_firstterm  double vector, real or complex
 *   niterations  nonnegative scalar integer value
 *   Nout         positive scalar integer value
 *
 * Output:
 *   Qfinal       complex double column vector of length Nout
 *
 * Equivalent MATLAB code:
 *
 *   Q = Q_firstterm;
 *   Qfinal = Q;
 *   for ii = 1:niterations
 *       Q = accumarray(out_idx, coeff .* Q(in_idx), [Nout 1]);
 *       Qfinal = Qfinal + Q;
 *   end
 *
 * Notes:
 *   - in_idx, out_idx, coeff must all have the same number of elements.
 *   - This uses the interleaved complex API.
 *   - Compile with:
 *         mex -R2018a EDinteq_iterate_mex.c
 */

static void require_vector(const mxArray *a, const char *name)
{
    if (!(mxGetM(a) == 1 || mxGetN(a) == 1)) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:notVector",
                          "%s must be a vector.", name);
    }
}

static mwSize get_vector_length(const mxArray *a, const char *name)
{
    require_vector(a, name);
    return mxGetNumberOfElements(a);
}

static uint64_t get_nonnegative_integer_scalar(const mxArray *a, const char *name)
{
    double v;
    uint64_t u;

    if (!mxIsNumeric(a) || mxIsComplex(a) || mxGetNumberOfElements(a) != 1) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:badScalar",
                          "%s must be a real numeric scalar.", name);
    }

    v = mxGetScalar(a);

    if (v < 0.0) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:negativeScalar",
                          "%s must be nonnegative.", name);
    }

    u = (uint64_t)v;
    if ((double)u != v) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:nonIntegerScalar",
                          "%s must be an integer-valued scalar.", name);
    }

    return u;
}

static mxComplexDouble get_coeff_at(mwSize k,
                                    bool coeff_is_complex,
                                    const mxComplexDouble *coeff_c,
                                    const double *coeff_r)
{
    mxComplexDouble a;
    if (coeff_is_complex) {
        a = coeff_c[k];
    } else {
        a.real = coeff_r[k];
        a.imag = 0.0;
    }
    return a;
}

static mxComplexDouble get_q0_at(mwSize k,
                                 bool q0_is_complex,
                                 const mxComplexDouble *q0_c,
                                 const double *q0_r)
{
    mxComplexDouble x;
    if (q0_is_complex) {
        x = q0_c[k];
    } else {
        x.real = q0_r[k];
        x.imag = 0.0;
    }
    return x;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    const mxArray *mx_in_idx;
    const mxArray *mx_out_idx;
    const mxArray *mx_coeff;
    const mxArray *mx_q0;

    const uint32_T *in_idx;
    const uint32_T *out_idx;

    const mxComplexDouble *coeff_c;
    const double *coeff_r;
    const mxComplexDouble *q0_c;
    const double *q0_r;

    bool coeff_is_complex;
    bool q0_is_complex;

    uint64_t niterations_u64;
    uint64_t Nout_u64;

    mwSize Nnz;
    mwSize Nq0;
    mwSize Nout;
    mwSize k;
    uint64_t iter;

    mxArray *out_mx;
    mxComplexDouble *Qfinal;
    mxComplexDouble *Qbuf1;
    mxComplexDouble *Qbuf2;
    mxComplexDouble *Qin;
    mxComplexDouble *Qout;
    mxComplexDouble *tmp;

    if (nrhs != 6) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:nrhs",
                          "Expected 6 inputs: in_idx, out_idx, coeff, Q_firstterm, niterations, Nout.");
    }

    if (nlhs > 1) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:nlhs",
                          "One output only.");
    }

    mx_in_idx  = prhs[0];
    mx_out_idx = prhs[1];
    mx_coeff   = prhs[2];
    mx_q0      = prhs[3];

    if (!mxIsUint32(mx_in_idx) || mxIsComplex(mx_in_idx)) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:badInIdx",
                          "in_idx must be a real uint32 vector.");
    }

    if (!mxIsUint32(mx_out_idx) || mxIsComplex(mx_out_idx)) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:badOutIdx",
                          "out_idx must be a real uint32 vector.");
    }

    if (!mxIsDouble(mx_coeff)) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:badCoeff",
                          "coeff must be a double vector, real or complex.");
    }

    if (!mxIsDouble(mx_q0)) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:badQ0",
                          "Q_firstterm must be a double vector, real or complex.");
    }

    Nnz = get_vector_length(mx_in_idx, "in_idx");

    if (get_vector_length(mx_out_idx, "out_idx") != Nnz) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:lengthMismatch",
                          "in_idx and out_idx must have the same number of elements.");
    }

    if (get_vector_length(mx_coeff, "coeff") != Nnz) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:lengthMismatch",
                          "coeff must have the same number of elements as in_idx.");
    }

    Nq0 = get_vector_length(mx_q0, "Q_firstterm");

    niterations_u64 = get_nonnegative_integer_scalar(prhs[4], "niterations");
    Nout_u64        = get_nonnegative_integer_scalar(prhs[5], "Nout");

    if (Nout_u64 == 0) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:badNout",
                          "Nout must be positive.");
    }

    if (Nout_u64 > (uint64_t)(~(mwSize)0)) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:tooLarge",
                          "Nout is too large for this platform.");
    }

    Nout = (mwSize)Nout_u64;

    if (Nq0 != Nout) {
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:sizeMismatch",
                          "numel(Q_firstterm) must equal Nout.");
    }

    in_idx  = (const uint32_T *)mxGetData(mx_in_idx);
    out_idx = (const uint32_T *)mxGetData(mx_out_idx);

    coeff_is_complex = mxIsComplex(mx_coeff);
    q0_is_complex    = mxIsComplex(mx_q0);

    if (coeff_is_complex) {
        coeff_c = mxGetComplexDoubles(mx_coeff);
        coeff_r = NULL;
        if (coeff_c == NULL) {
            mexErrMsgIdAndTxt("EDinteq_iterate_mex:dataAccess",
                              "Failed to access complex coeff data.");
        }
    } else {
        coeff_r = mxGetDoubles(mx_coeff);
        coeff_c = NULL;
        if (coeff_r == NULL) {
            mexErrMsgIdAndTxt("EDinteq_iterate_mex:dataAccess",
                              "Failed to access real coeff data.");
        }
    }

    if (q0_is_complex) {
        q0_c = mxGetComplexDoubles(mx_q0);
        q0_r = NULL;
        if (q0_c == NULL) {
            mexErrMsgIdAndTxt("EDinteq_iterate_mex:dataAccess",
                              "Failed to access complex Q_firstterm data.");
        }
    } else {
        q0_r = mxGetDoubles(mx_q0);
        q0_c = NULL;
        if (q0_r == NULL) {
            mexErrMsgIdAndTxt("EDinteq_iterate_mex:dataAccess",
                              "Failed to access real Q_firstterm data.");
        }
    }

    for (k = 0; k < Nnz; ++k) {
        if (in_idx[k] < 1 || in_idx[k] > (uint32_T)Nout) {
            mexErrMsgIdAndTxt("EDinteq_iterate_mex:inIdxRange",
                              "in_idx contains values outside 1..Nout.");
        }

        if (out_idx[k] < 1 || out_idx[k] > (uint32_T)Nout) {
            mexErrMsgIdAndTxt("EDinteq_iterate_mex:outIdxRange",
                              "out_idx contains values outside 1..Nout.");
        }
    }

    out_mx = mxCreateDoubleMatrix(Nout, 1, mxCOMPLEX);
    plhs[0] = out_mx;
    Qfinal = mxGetComplexDoubles(out_mx);

    Qbuf1 = (mxComplexDouble *)mxCalloc(Nout, sizeof(mxComplexDouble));
    Qbuf2 = (mxComplexDouble *)mxCalloc(Nout, sizeof(mxComplexDouble));

    if (Qbuf1 == NULL || Qbuf2 == NULL) {
        if (Qbuf1 != NULL) mxFree(Qbuf1);
        if (Qbuf2 != NULL) mxFree(Qbuf2);
        mexErrMsgIdAndTxt("EDinteq_iterate_mex:alloc",
                          "Failed to allocate work buffers.");
    }

    /* Initialize Q = Q_firstterm and Qfinal = Q_firstterm */
    for (k = 0; k < Nout; ++k) {
        mxComplexDouble x = get_q0_at(k, q0_is_complex, q0_c, q0_r);
        Qbuf1[k] = x;
        Qfinal[k] = x;
    }

    Qin  = Qbuf1;
    Qout = Qbuf2;

    for (iter = 0; iter < niterations_u64; ++iter) {
        memset(Qout, 0, Nout * sizeof(mxComplexDouble));

        for (k = 0; k < Nnz; ++k) {
            mwSize i_in  = (mwSize)(in_idx[k]  - 1);
            mwSize i_out = (mwSize)(out_idx[k] - 1);

            mxComplexDouble a = get_coeff_at(k, coeff_is_complex, coeff_c, coeff_r);
            mxComplexDouble x = Qin[i_in];

            Qout[i_out].real += a.real * x.real - a.imag * x.imag;
            Qout[i_out].imag += a.real * x.imag + a.imag * x.real;
        }

        for (k = 0; k < Nout; ++k) {
            Qfinal[k].real += Qout[k].real;
            Qfinal[k].imag += Qout[k].imag;
        }

        tmp  = Qin;
        Qin  = Qout;
        Qout = tmp;
    }

    mxFree(Qbuf1);
    mxFree(Qbuf2);
}
