/*=================================================================
 *
 * ESIE2calcedgeinteqmatrixsub2_mex.c	
 *
 *	Peter Svensson 12 February 2013
 *
 * The calling syntax is:
 *
 * [ivuse,n1hormat,n2vertmat,n3vertmat] = ESIE2calcedgeinteqmatrixsub2_mex(n1,n2,3);
 *
 *=================================================================*
 * First version 12 February 2013
 * 1 April 2015 Changed name to ESIE2calcedgeinteqmatrixsub2_mex
 */
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	N1	prhs[0]    /* Integer, uint16 */
#define	N2	prhs[1]    /* Integer, uint16 */
#define	N3	prhs[2]    /* Integer, uint16 */

/* Output Arguments */

#define	IVUSE               plhs[0]     /* Array of integers, uint32 */
#define	N1HORMAT			plhs[1]     /* Array of integers, uint32 */
#define	N2VERTMAT			plhs[2]     /* Array of integers, uint32 */
#define	N3VERTMAT			plhs[3]     /* Array of integers, uint32 */ 

/*******************************************************************************************
*******************************************************************************************
*******************************************************************************************
*******************************************************************************************
*******************************************************************************************/

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    
    double *ivuse,*n1hormat, *n2vertmat, *n3vertmat;
    double ivusevalue;
/*    double *qoutput_real, *qoutput_imag, *pr;
    double multres_real,multres_imag;  */
    
    unsigned int n1, n2, n3, ntot;
    unsigned int i,j,k,counter,counter2;
        
	/********************************************/
	/*                                          */
    /* Check for proper number of arguments     */
	/*                                          */
    
    if (nrhs != 3) { 
        printf("   nrhs = %d\n",nrhs);
		mexErrMsgTxt("Three input arguments required: n1, n2, n3"); 
    } else if (nlhs != 4) {
        printf("   nlhs = %d\n",nlhs);
		mexErrMsgTxt("Four output arguments required: ivuse,n1hormat,n2vertmat,n3vertmat");
    } 
        
	/********************************************/
	/*                                          */
	/* Input parameters                         */
	/*                                          */
    
    /* Get the input value  */

    n1 = (int)(mxGetScalar(N1));
    n2 = (int)(mxGetScalar(N2));
    n3 = (int)(mxGetScalar(N3));
    
    ntot = n1*n2*n3;
                  
	/********************************************/
	/*                                          */
	/* Output parameters                        */
	/* Qoutput is assigned further down, when		*/
	/* we know which length it should be        */

    
/***    QOUTPUT = mxCreateDoubleMatrix(nqvec, 1, mxCOMPLEX); 
	qoutput_real = mxGetPr(QOUTPUT);
	qoutput_imag = mxGetPi(QOUTPUT);
/*    ivuse = mxGetPr(IVUSE);
    n1hormat = mxGetPr(IVUSE);
    n2vertmat = mxGetPr(IVUSE);
    n3vertmat = mxGetPr(IVUSE);*/
    
    IVUSE = mxCreateDoubleMatrix(ntot, 1, mxREAL); 
	ivuse = mxGetPr(IVUSE);
    N1HORMAT = mxCreateDoubleMatrix(ntot, 1, mxREAL); 
	n1hormat = mxGetPr(N1HORMAT);
    N2VERTMAT = mxCreateDoubleMatrix(ntot, 1, mxREAL); 
	n2vertmat = mxGetPr(N2VERTMAT);
    N3VERTMAT = mxCreateDoubleMatrix(ntot, 1, mxREAL); 
	n3vertmat = mxGetPr(N3VERTMAT);
    

	/********************************************/
	/*                                          */
	/* The actual calculations                  */
	/*                                          */

/*    printf("ntotal: %d\n", ntot);
    
/*    printf("qinput value 1 = %f +j* %f\n",qinput_real[0],qinput_imag[0]); */
    
    /* Fill the output vector with the first contribution: the input vector*/
    /* Also create a temporary q-vector for storage inside the for-loop */
    for (i=0;i<ntot;++i) {
		ivuse[i] = 0;
        n1hormat[i] = i;
        n2vertmat[i] = i+1;
        n3vertmat[i] = i*i;
	}
    
    counter = 0;
    counter2 = 0;
    ivusevalue = 1;
    for (i=0;i<n2;++i) {
        for (j=0;j<n1;++j) {
            for (k=0;k<n3;++k) {
                n1hormat[counter] = j+1;
                n2vertmat[counter] = i+1;
                n3vertmat[counter] = k+1;
                ivuse[counter] = ivusevalue;
                counter = counter +1;
                counter2 = counter2 + 1;
                if (counter2 < n1*n3) {
                    ivusevalue = ivusevalue + n2;
                } else {
                    ivusevalue = ivusevalue + n2+1;
                    counter2 = 0;
                }
            }
        }
        
        
    }
    
    
    
    
    

/*    printf("value list size: %d\n",nvaluevec);

        printf("   No. iteration: %d\n",niterations);*/

/*    printf("inputlocvec pos. 1 = %d\n",inputlocvector[0]);***/
    
    
/*    for (counter=0;counter<niterations;++counter) {


        for (i=0;i<nvaluevec;++i) {
            multres_real = valuevector_real[i]*qackum_real[inputlocvector[i]] - valuevector_imag[i]*qackum_imag[inputlocvector[i]];
            multres_imag = valuevector_real[i]*qackum_imag[inputlocvector[i]] + valuevector_imag[i]*qackum_real[inputlocvector[i]];
            
            qtemp_real[outputlocvector[i]] = qtemp_real[outputlocvector[i]] + multres_real;
            qtemp_imag[outputlocvector[i]] = qtemp_imag[outputlocvector[i]] + multres_imag;
                    
        }
        for (i=0;i<nqvec;++i) {
        	qoutput_real[i] = qoutput_real[i] + qtemp_real[i];
        	qoutput_imag[i] = qoutput_imag[i] + qtemp_imag[i];
            qackum_real[i] = qtemp_real[i];
            qackum_imag[i] = qtemp_imag[i];
            qtemp_real[i] = 0;
            qtemp_imag[i] = 0;
    	}
        
    }*/ 
    	
    return;
    
}
