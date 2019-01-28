/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2018 The MathWorks, Inc.
 *
 *========================================================*/

#include "mex.h"
#include "string.h"

mwSize U_cor(mwSize i, mwSize j, mwSize N){
    mwSize U;
    U = i-1+(j-1)*(N-1);
    return U;
}

mwSize V_cor(mwSize i, mwSize j, mwSize N){
    mwSize V;
    V = j-1+(i-1)*(N-1);
    return V;
}

mwSize P_cor(mwSize i, mwSize j, mwSize N){
    mwSize P;
    P = i-1+(j-1)*N;
    return P;
}

/* The computational routine */
void testMex(mwSize N, double *inU, double *inV, double *inP, double *outU, double *outV,double *outP, int nrow1, int nrow2, int nrow3)
{
    mwSize i,j;
    double r,delta,h;
    h = 1./((double)N);
    /* copy parameter */
    memcpy(outU,inU,sizeof(double)*nrow1);
    memcpy(outV,inV,sizeof(double)*nrow2);
    memcpy(outP,inP,sizeof(double)*nrow3);

    for (i=2; i<N; i++){
        for (j=2; j<N; j++){
            r = (outU[U_cor(i,j,N)]-outU[U_cor(i-1,j,N)]+outV[V_cor(i,j,N)]-outV[V_cor(i,j-1,N)])/h;
            delta = r*h/4.;
            outU[U_cor(i-1,j,N)] = outU[U_cor(i-1,j,N)]+delta;
            outU[U_cor(i,j,N)] = outU[U_cor(i,j,N)]-delta;
            outV[V_cor(i,j-1,N)] = outV[V_cor(i,j-1,N)]+delta;
            outV[V_cor(i,j,N)] = outV[V_cor(i,j,N)]-delta;
            outP[P_cor(i,j,N)] = outP[P_cor(i,j,N)]-r;
            outP[P_cor(i+1,j,N)] = outP[P_cor(i+1,j,N)]+r/4.;
            outP[P_cor(i-1,j,N)] = outP[P_cor(i-1,j,N)]+r/4.;
            outP[P_cor(i,j+1,N)] = outP[P_cor(i,j+1,N)]+r/4.;
            outP[P_cor(i,j-1,N)] = outP[P_cor(i,j-1,N)]+r/4.;
        }
    }

    i = 1;
    for (j=2; j<N; j++){
        r = (outU[U_cor(i,j,N)]+outV[V_cor(i,j,N)]-outV[V_cor(i,j-1,N)])/h;
        delta = r*h/3.;
        outU[U_cor(i,j,N)] = outU[U_cor(i,j,N)]-delta;
        outV[V_cor(i,j-1,N)] = outV[V_cor(i,j-1,N)]+delta;
        outV[V_cor(i,j,N)] = outV[V_cor(i,j,N)]-delta;
        outP[P_cor(i,j,N)] = outP[P_cor(i,j,N)]-4.*r/3.;
        outP[P_cor(i+1,j,N)] = outP[P_cor(i+1,j,N)]+r/3.;
        outP[P_cor(i,j+1,N)] = outP[P_cor(i,j+1,N)]+r/3.;
        outP[P_cor(i,j-1,N)] = outP[P_cor(i,j-1,N)]+r/3.;
    }

    i = N;
    for (j=2; j<N; j++){
        r = (-outU[U_cor(i-1,j,N)]+outV[V_cor(i,j,N)]-outV[V_cor(i,j-1,N)])/h;
        delta = r*h/3.;
        outU[U_cor(i-1,j,N)] = outU[U_cor(i-1,j,N)]+delta;
        outV[V_cor(i,j-1,N)] = outV[V_cor(i,j-1,N)]+delta;
        outV[V_cor(i,j,N)] = outV[V_cor(i,j,N)]-delta;
        outP[P_cor(i,j,N)] = outP[P_cor(i,j,N)]-4.*r/3.;
        outP[P_cor(i-1,j,N)] = outP[P_cor(i-1,j,N)]+r/3.;
        outP[P_cor(i,j+1,N)] = outP[P_cor(i,j+1,N)]+r/3.;
        outP[P_cor(i,j-1,N)] = outP[P_cor(i,j-1,N)]+r/3.;
    }

    j = 1;
    for (i=2; i<N; i++){
        r = (outU[U_cor(i,j,N)]-outU[U_cor(i-1,j,N)]+outV[V_cor(i,j,N)])/h;
        delta = r*h/3.;
        outU[U_cor(i-1,j,N)] = outU[U_cor(i-1,j,N)]+delta;
        outU[U_cor(i,j,N)] = outU[U_cor(i,j,N)]-delta;
        outV[V_cor(i,j,N)] = outV[V_cor(i,j,N)]-delta;
        outP[P_cor(i,j,N)] = outP[P_cor(i,j,N)]-4.*r/3.;
        outP[P_cor(i+1,j,N)] = outP[P_cor(i+1,j,N)]+r/3.;
        outP[P_cor(i-1,j,N)] = outP[P_cor(i-1,j,N)]+r/3.;
        outP[P_cor(i,j+1,N)] = outP[P_cor(i,j+1,N)]+r/3.;
    }

    j = N;
    for (i=2; i<N; i++){
        r = (outU[U_cor(i,j,N)]-outU[U_cor(i-1,j,N)]-outV[V_cor(i,j-1,N)])/h;
        delta = r*h/3;
        outU[U_cor(i-1,j,N)] = outU[U_cor(i-1,j,N)]+delta;
        outU[U_cor(i,j,N)] = outU[U_cor(i,j,N)]-delta;
        outV[V_cor(i,j-1,N)] = outV[V_cor(i,j-1,N)]+delta;
        outP[P_cor(i,j,N)] = outP[P_cor(i,j,N)]-4.*r/3.;
        outP[P_cor(i+1,j,N)] = outP[P_cor(i+1,j,N)]+r/3.;
        outP[P_cor(i-1,j,N)] = outP[P_cor(i-1,j,N)]+r/3.;
        outP[P_cor(i,j-1,N)] = outP[P_cor(i,j-1,N)]+r/3.;
    }
    
    i = 1; j = 1;
    r = (outU[U_cor(i,j,N)]+outV[V_cor(i,j,N)])/h;
    delta = r*h/2;
    outU[U_cor(i,j,N)] = outU[U_cor(i,j,N)]-delta;
    outV[V_cor(i,j,N)] = outV[V_cor(i,j,N)]-delta;
    outP[P_cor(i,j,N)] = outP[P_cor(i,j,N)]-2.*r;
    outP[P_cor(i+1,j,N)] = outP[P_cor(i+1,j,N)]+r/2.;
    outP[P_cor(i,j+1,N)] = outP[P_cor(i,j+1,N)]+r/2.;

    i = 1; j = N;
    r = (outU[U_cor(i,j,N)]-outV[V_cor(i,j-1,N)])/h;
    delta = r*h/2;
    outU[U_cor(i,j,N)] = outU[U_cor(i,j,N)]-delta;
    outV[V_cor(i,j-1,N)] = outV[V_cor(i,j-1,N)]+delta;
    outP[P_cor(i,j,N)] = outP[P_cor(i,j,N)]-2.*r;
    outP[P_cor(i+1,j,N)] = outP[P_cor(i+1,j,N)]+r/2.;
    outP[P_cor(i,j-1,N)] = outP[P_cor(i,j-1,N)]+r/2.;

    i = N; j = 1;
    r = (-outU[U_cor(i-1,j,N)]+outV[V_cor(i,j,N)])/h;
    delta = r*h/2;
    outU[U_cor(i-1,j,N)] = outU[U_cor(i-1,j,N)]+delta;
    outV[V_cor(i,j,N)] = outV[V_cor(i,j,N)]-delta;
    outP[P_cor(i,j,N)] = outP[P_cor(i,j,N)]-2.*r;
    outP[P_cor(i-1,j,N)] = outP[P_cor(i-1,j,N)]+r/2.;
    outP[P_cor(i,j+1,N)] = outP[P_cor(i,j+1,N)]+r/2.;

    i = N; j = N;
    r = (-outU[U_cor(i-1,j,N)]-outV[V_cor(i,j-1,N)])/h;
    delta = r*h/2.;
    outU[U_cor(i-1,j,N)] = outU[U_cor(i-1,j,N)]+delta;
    outV[V_cor(i,j-1,N)] = outV[V_cor(i,j-1,N)]+delta;
    outP[P_cor(i,j,N)] = outP[P_cor(i,j,N)]-2.*r;
    outP[P_cor(i-1,j,N)] = outP[P_cor(i-1,j,N)]+r/2.;
    outP[P_cor(i,j-1,N)] = outP[P_cor(i,j-1,N)]+r/2.;
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    size_t nrow1, nrow2, nrow3;
    int N;
    double *inU, *inV, *inP;               /* 1xN^2 input matrix */
    double *outU, *outV, *outP;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Four inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Three output required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
    }
    
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    
    /* check that number of rows in second input argument is 1 */
    if(mxGetN(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notColVector","Input must be a column vector.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetN(prhs[2])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notColVector","Input must be a column vector.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetN(prhs[3])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notColVector","Input must be a column vector.");
    }
    
    /* get the value of the scalar input  */
    N = mxGetScalar(prhs[0]);

    nrow1 = mxGetM(prhs[1]);
    nrow2 = mxGetM(prhs[2]);
    nrow3 = mxGetM(prhs[3]);

    if( (int)nrow1!= N*(N-1) || 
        (int)nrow2!= N*(N-1) ||
        (int)nrow3!= N*N) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Wrong Size");
    }

    /* create a pointer to the real data in the input matrix  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    inU = mxGetDoubles(prhs[1]);
    inV = mxGetDoubles(prhs[2]);
    inP = mxGetDoubles(prhs[3]);
    #else
    inU = mxGetPr(prhs[1]);
    inV = mxGetPr(prhs[2]);
    inP = mxGetPr(prhs[3]);
    #endif


    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nrow1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)nrow2,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)nrow3,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    outU = mxGetDoubles(plhs[0]);
    outV = mxGetDoubles(plhs[1]);
    outP = mxGetDoubles(plhs[2]);
    #else
    outU = mxGetPr(plhs[0]);
    outV = mxGetPr(plhs[1]);
    outP = mxGetPr(plhs[2]);
    #endif

    /* call the computational routine */
    testMex((mwSize)N, inU, inV, inP, outU, outV, outP,(int)nrow1, (int) nrow2, (int)nrow3);
}
