#include "GPUComplex.h"

/*
 * Return the square of a complex number 
 */
__device__ inline GPUComplex d_GPUComplex_square(GPUComplex& src) {
    return GPUComplex(src.x*src.x - src.y*src.y, 2*src.x*src.y);
}

/*
 * Return the conjugate of a complex number 
 */
__device__ inline GPUComplex d_GPUComplex_conj(const GPUComplex& src) {
return GPUComplex(src.x, -src.y);
}


/*
 * Return the product of 2 complex numbers 
 */
__device__ inline GPUComplex d_GPUComplex_product(const GPUComplex& a, const GPUComplex& b) {
    return GPUComplex(a.x * b.x - a.y*b.y, a.x * b.y + a.y*b.x);
}


/*
 * Return the absolute of a complex number 
 */
__device__ inline double d_GPUComplex_abs(const GPUComplex& src) {
    return sqrt(src.x * src.x + src.y * src.y);
}

/*
 *  result = a * b * c (a = complex ; b,c = double) 
 */
__device__ inline GPUComplex d_GPUComplex_mult(GPUComplex& a, double b, double c) {
    return GPUComplex(a.x * b * c, a.y * b * c);
}

/*
 * Return the complex number c = a * b (a is complex, b is double) 
 */
__device__ inline GPUComplex d_GPUComplex_mult(const GPUComplex& a, double b) {
   return GPUComplex(a.x*b, a.y*b);

}

/*
 * Return the complex number a += b * c  
 */
__device__ inline void d_GPUComplex_fma(GPUComplex& a, const GPUComplex& b, const GPUComplex& c) {
    a.x += b.x * c.x - b.y*c.y ;
    a.y += b.x * c.y + b.y*c.x ;
}

/*
 * Return the complex number a -= b * c  
 */
__device__ inline void d_GPUComplex_fms(GPUComplex& a, const GPUComplex& b, const GPUComplex& c) {
    a.x -= b.x * c.x - b.y*c.y ;
    a.y -= b.x * c.y + b.y*c.x ;
}


__device__ inline GPUComplex d_doubleMinusGPUComplex(const double &a, GPUComplex& src) {
    return GPUComplex(a-src.x, -src.y);
}

__device__ inline GPUComplex d_doublePlusGPUComplex(double a, GPUComplex& src) {
    return GPUComplex(a+src.x, src.y);
}

__device__ inline double d_GPUComplex_real( const GPUComplex& src) {
    return src.x;
}

__device__ inline double d_GPUComplex_imag( const GPUComplex& src) {
    return src.y;
}

__device__ inline void d_GPUComplex_plusEquals( GPUComplex& a, const GPUComplex & b) {
    a.x += b.x;
    a.y += b.y;
}

__device__ void inline d_GPUComplex_Equals( GPUComplex& a, const GPUComplex & b) {
    a.x = b.x;
    a.y = b.y;
}

__device__ void d_print( const GPUComplex& a) {
    printf("( %f, %f) ", a.x, a.y);
    printf("\n");
}

__device__ void ncoulsKernel(GPUComplex& mygpvar1, GPUComplex& wdiff, GPUComplex& aqsntemp_index, GPUComplex& wtilde_array_index, GPUComplex& I_eps_array_index, double vcoul_igp, double& achtemp_re_loc, double& achtemp_im_loc)
{
    double rden = 1/(wdiff.x*wdiff.x + wdiff.y*wdiff.y);

    achtemp_re_loc += d_GPUComplex_product(d_GPUComplex_product(mygpvar1, aqsntemp_index),\
        d_GPUComplex_product(d_GPUComplex_mult(d_GPUComplex_product(wtilde_array_index, d_GPUComplex_conj(wdiff)), rden), I_eps_array_index)).x * 0.5 * vcoul_igp;
    achtemp_im_loc += d_GPUComplex_product(d_GPUComplex_product(mygpvar1, aqsntemp_index),\
        d_GPUComplex_product(d_GPUComplex_mult(d_GPUComplex_product(wtilde_array_index, d_GPUComplex_conj(wdiff)), rden), I_eps_array_index)).y * 0.5 * vcoul_igp;

}

__global__ void d_flagOCC_solver(double *wx_array, GPUComplex *wtilde_array, GPUComplex* asxtemp, GPUComplex *aqsmtemp, GPUComplex *aqsntemp, GPUComplex *I_eps_array, int* inv_igp_index, int* indinv, int ncouls, int nvband, int ngpown, double* vcoul)
{
    int n1 = blockIdx.x ;
    int my_igp = blockIdx.y;

    if(n1 < nvband && my_igp < ngpown)
    {
        int loopOverncouls = 1, \
            threadsPerBlock = 128;

        if(ncouls > threadsPerBlock)
            loopOverncouls = ncouls / threadsPerBlock;

        for(int iw = nstart; iw < nend; ++iw)
        {
            double wxt = wx_array[iw];
            for( int x = 0; x < loopOverncouls && threadIdx.x < threadsPerBlock ; ++x)
            {
                int indigp = inv_igp_index[my_igp];
                int igp = indinv[indigp];
                GPUComplex ssxt = GPUComplex(0.00, 0.00);
                GPUComplex scht = GPUComplex(0.00, 0.00);
                        
                GPUComplex expr = GPUComplex(0.50, 0.50);
                GPUComplex expr0 = GPUComplex(0.00, 0.00);
                GPUComplex matngmatmgp = GPUComplex(0.00, 0.00);
                GPUComplex matngpmatmg = GPUComplex(0.00, 0.00);
                GPUComplex wtilde(0.00, 0.00);
                GPUComplex wtilde2(0.00, 0.00);
                GPUComplex Omega2(0.00, 0.00);
                GPUComplex mygpvar1(0.00, 0.00);
                GPUComplex mygpvar2(0.00, 0.00);
                GPUComplex ssx(0.00, 0.00);
        
                for(int ig=0; ig<ncouls; ++ig)
                {
                    d_GPUComplex_Equals( wtilde , wtilde_array[my_igp*ncouls+ig]);
                    d_GPUComplex_Equals( wtilde2 , d_GPUComplex_square(wtilde));
                    d_GPUComplex_Equals( Omega2 , d_GPUComplex_product(wtilde2,I_eps_array[my_igp*ncouls+ig]));
                    d_GPUComplex_Equals( mygpvar1 , d_GPUComplex_conj(aqsmtemp[n1*ncouls+igp]));
                    d_GPUComplex_Equals( mygpvar2 , aqsmtemp[n1*ncouls+igp]);
                    d_GPUComplex_Equals( matngmatmgp , d_GPUComplex_product(aqsntemp[n1*ncouls+ig] , mygpvar1));
                    if(ig != igp) d_GPUComplex_Equals( matngpmatmg , d_GPUComplex_product(d_GPUComplex_conj(aqsmtemp[n1*ncouls+ig]) , mygpvar2));
        
                    double to1 = 1e-6;
                    double sexcut = 4.0;
                    double limitone = 1.0/(to1*4.0);
                    double limittwo = pow(0.5,2);
                
                    GPUComplex wdiff = d_doubleMinusGPUComplex(wxt , wtilde);
                
                    double rden = 1/d_GPUComplex_real(d_GPUComplex_product(d_doubleMinusGPUComplex(wxt , wtilde), d_GPUComplex_conj(d_doubleMinusGPUComplex(wxt , wtilde))));
                    GPUComplex delw = d_GPUComplex_mult(d_GPUComplex_product(wtilde , d_GPUComplex_conj(d_doubleMinusGPUComplex(wxt , wtilde))) , rden);
                    double delwr = d_GPUComplex_real(d_GPUComplex_product(delw , d_GPUComplex_conj(delw)));
                    double wdiffr = d_GPUComplex_real(d_GPUComplex_product(d_doubleMinusGPUComplex(wxt , wtilde), d_GPUComplex_conj(d_doubleMinusGPUComplex(wxt , wtilde))));
                
                    if((wdiffr > limittwo) && (delwr < limitone))
                    {
                       double cden = std::pow(wxt,2);
                        rden = std::pow(cden,2);
                        rden = 1.00 / rden;
                        d_GPUComplex_Equals(ssx , d_GPUComplex_mult(Omega2 , cden , rden));
                    }
                    else if (delwr > to1)
                    {
                        GPUComplex cden(0.00, 0.00);
                        d_GPUComplex_Equals(cden , d_GPUComplex_mult(d_GPUComplex_product(wtilde2, d_doublePlusGPUComplex((double)0.50, delw)), 4.00));
                       rden = d_GPUComplex_real(d_GPUComplex_product(cden , d_GPUComplex_conj(cden)));
                        rden = 1.00/rden;
                        d_GPUComplex_Equals(ssx , d_GPUComplex_product(d_GPUComplex_product(GPUComplex(-Omega2.x, -Omega2.y) , d_GPUComplex_conj(cden)), d_GPUComplex_mult(delw, rden)));
                    }
                    else
                    {
                        d_GPUComplex_Equals(ssx , expr0);
                    }
                
                    double ssxcutoff = d_GPUComplex_abs(I_eps_array[my_igp*ngpown+ig]) * sexcut;
                    if((d_GPUComplex_abs(ssx) > ssxcutoff) && (wxt < 0.00)) d_GPUComplex_Equals(ssx , expr0);
        
                    d_GPUComplex_plusEquals(ssxt, d_GPUComplex_product(matngmatmgp , ssxt));
                    d_GPUComplex_plusEquals(scht, d_GPUComplex_product(matngmatmgp , scht));
                }
                        d_GPUComplex_plusEquals(asxtemp[iw] , d_GPUComplex_mult(ssxt , vcoul[igp]));
            }
        }
    }
}

__global__  void NumberBands_kernel( GPUComplex *wtilde_array, GPUComplex *aqsntemp, GPUComplex* aqsmtemp, GPUComplex *I_eps_array, int ncouls, int ngpown, int number_bands, double* wx_array, double* achtemp_re, double* achtemp_im, double* vcoul, int* indinv, int* inv_igp_index, int numThreadsPerBlock, int stride)
{
    int n1 = blockIdx.x ;

    if(n1 < number_bands)
    {
        int loopOverngpown = 1, leftOverngpown = 0, \
            loopCounter = numThreadsPerBlock;

        if(ngpown > loopCounter)
        {
            loopOverngpown = ngpown / loopCounter;
            leftOverngpown = ngpown % loopCounter;
        }

        double achtemp_re_loc[nend-nstart];
        double achtemp_im_loc[nend-nstart];
        for(int iw = nstart; iw < nend; ++iw) {achtemp_re_loc[iw] = 0.00; achtemp_im_loc[iw] = 0.00;}

        for( int x = 0; x < loopOverngpown && threadIdx.x < loopCounter; ++x)
        {
            int my_igp = x*loopCounter + threadIdx.x;
        
            if(my_igp < ngpown)
            {
                int indigp = inv_igp_index[my_igp];
                int igp = indinv[indigp];

                //for(int iw = nstart; iw < nend; ++iw)
                //{
                    if(stride == 0)
                    {
                        for(int ig = 0; ig < ncouls; ++ig) 
                        { 
                            for(int iw = nstart; iw < nend; ++iw)
			    {
                            GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls+igp]);
                            GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                            ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
			    }
                        } //ncouls
                    }
                    else
                    {
                        for(int igmin = 0; igmin < stride; ++igmin)
                        {
                            for(int ig = 0; ig < ncouls; ig+=stride) 
                            {
                                for(int iw = nstart; iw < nend; ++iw)
			        {
                                GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls+igp]);
                                GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                                ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
			        }
                            }
                        }
                    } //else
                //} // iw
            } // ngpown
        }

        if(leftOverngpown)
        {
            int my_igp = loopOverngpown*loopCounter + threadIdx.x;
            if(my_igp < ngpown)
            {
                int indigp = inv_igp_index[my_igp];
                int igp = indinv[indigp];


                    for(int ig = 0; ig < ncouls; ++ig) 
                    { 
                for(int iw = nstart; iw < nend; ++iw)
                {
                        GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls+igp]);
                        GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                        ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);

                    } //ncouls

                } // iw
            } // ngpown
        }

        for(int iw = nstart; iw < nend; ++iw)
        {
            atomicAdd(&achtemp_re[iw] , achtemp_re_loc[iw] );
            atomicAdd(&achtemp_im[iw] , achtemp_im_loc[iw] );
        }
    }
}

__global__  void cudaNgpown_kernel( int n1, GPUComplex *wtilde_array, GPUComplex *aqsntemp, GPUComplex* aqsmtemp, GPUComplex *I_eps_array, int ncouls, int ngpown, int number_bands, double* wx_array, double* achtemp_re, double* achtemp_im, double* vcoul, int* indinv, int* inv_igp_index, int numThreadsPerBlock)
{
    int my_igp = blockIdx.x;

    if(my_igp < ngpown )
    {
        int loopOverncouls = 1, leftOverncouls = 0;

        if(ncouls > numThreadsPerBlock)
        {
            loopOverncouls = ncouls / numThreadsPerBlock;
            leftOverncouls = ncouls % numThreadsPerBlock;
        }

        double achtemp_re_loc[nend-nstart], achtemp_im_loc[nend-nstart];
        int indigp = inv_igp_index[my_igp];
        int igp = indinv[indigp];

        for(int iw = nstart; iw < nend; ++iw)
        {
            achtemp_re_loc[iw] = 0.00; achtemp_im_loc[iw] = 0.00;
            for( int x = 0; x < loopOverncouls && threadIdx.x < numThreadsPerBlock ; ++x)
            { 
                int ig = x*numThreadsPerBlock + threadIdx.x;
                if(ig < ncouls)
                {
                    GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls +igp]);
                    GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                    ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
                }
            }
            if(leftOverncouls)
            {
                int ig = loopOverncouls*numThreadsPerBlock + threadIdx.x;
                if(ig < ncouls)
                {
                    GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls +igp]);
                    GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                    ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
                }
            }

            atomicAdd(&achtemp_re[iw] , achtemp_re_loc[iw] );
            atomicAdd(&achtemp_im[iw] , achtemp_im_loc[iw] );
        }
    }
}

__global__  void NumBandNgpown_kernel( GPUComplex *wtilde_array, GPUComplex *aqsntemp, GPUComplex* aqsmtemp, GPUComplex *I_eps_array, int ncouls, int ngpown, int number_bands, double* wx_array, double* achtemp_re, double* achtemp_im, double* vcoul, int* indinv, int* inv_igp_index, int numThreadsPerBlock, int stride)
{
    int n1 = blockIdx.x;
    int my_igp = blockIdx.y;

    if((n1 < number_bands ) && (my_igp < ngpown) )
    {
        int loopOverncouls = 1, leftOverncouls = 0;
        if(ncouls > numThreadsPerBlock)
        {
            loopOverncouls = ncouls / numThreadsPerBlock;
            leftOverncouls = ncouls % numThreadsPerBlock;
        }

        double achtemp_re_loc[nend-nstart], achtemp_im_loc[nend-nstart];
        int indigp = inv_igp_index[my_igp];
        int igp = indinv[indigp];
        for(int iw = nstart; iw < nend; ++iw) {achtemp_re_loc[iw] = 0.00; achtemp_im_loc[iw] = 0.00;}

        if(stride == 0)
        {
            for( int x = 0; x < loopOverncouls && threadIdx.x < numThreadsPerBlock ; ++x)
            { 
                int ig = x*numThreadsPerBlock + threadIdx.x;
                if(ig < ncouls)
                {
                    for(int iw = nstart; iw < nend; ++iw)
                    {
                        GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls +igp]);
                        GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                        ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
                    }
                }
            }
            if(leftOverncouls)
            {
                int ig = loopOverncouls*numThreadsPerBlock + threadIdx.x;
                if(ig < ncouls)
                {
                    for(int iw = nstart; iw < nend; ++iw)
                    {
                        GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls +igp]);
                        GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                        ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
                   }
                }
            }
        }
        else
        {
            for(int igmin = 0; igmin < stride; ++igmin)
            { 
                for( int x = 0; x < loopOverncouls/stride && threadIdx.x < numThreadsPerBlock ; ++x)
                {
                    int ig = (x*numThreadsPerBlock + threadIdx.x) * stride + igmin ;
                    if(ig < ncouls)
                    {
                      for(int iw = nstart; iw < nend; ++iw)
                      {
                          GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls +igp]);
                          GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                          ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
                      }
                    }
                }
            }
            if(leftOverncouls)
            {
                for(int igmin = 0; igmin < stride; ++igmin)
                {
                    int ig = loopOverncouls*numThreadsPerBlock + threadIdx.x*stride + igmin;
                    if(ig < ncouls)
                    {
                      for(int iw = nstart; iw < nend; ++iw)
                      {
                          GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls +igp]);
                          GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                          ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
                      }
                    }
                }
            }
        }

        for(int iw = nstart; iw < nend; ++iw)
        {
            atomicAdd(&achtemp_re[iw] , achtemp_re_loc[iw] );
            atomicAdd(&achtemp_im[iw] , achtemp_im_loc[iw] );
        }
    }
}

__global__  void NgpownNcouls_kernel( GPUComplex *wtilde_array, GPUComplex *aqsntemp, GPUComplex* aqsmtemp, GPUComplex *I_eps_array, int ncouls, int ngpown, int number_bands, double* wx_array, double* achtemp_re, double* achtemp_im, double* vcoul, int* indinv, int* inv_igp_index, int numThreadsPerBlock)
{
    int n1 = blockIdx.x ;

    if(n1 < number_bands)
    {
        int loopOverngpown = 1, leftOverngpown = 0, \
            loopCounter = numThreadsPerBlock;

        if(ngpown > loopCounter)
        {
            loopOverngpown = ngpown / loopCounter;
            leftOverngpown = ngpown % loopCounter;
        }

        double achtemp_re_loc[nend-nstart];
        double achtemp_im_loc[nend-nstart];
        for(int iw = nstart; iw < nend; ++iw) {achtemp_re_loc[iw] = 0.00; achtemp_im_loc[iw] = 0.00;}

        for( int x = 0; x < loopOverngpown && threadIdx.x < loopCounter; ++x)
        {
            int my_igp = x*loopCounter + threadIdx.x;
        
            if(my_igp < ngpown)
            {
                int indigp = inv_igp_index[my_igp];
                int igp = indinv[indigp];

                for(int iw = nstart; iw < nend; ++iw)
                {

                    for(int ig = 0; ig < ncouls; ++ig) 
                    { 
                        GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls+igp]);
                        GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                        ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
                    } //ncouls
                } // iw
            } // ngpown
        }

        if(leftOverngpown)
        {
            int my_igp = loopOverngpown*loopCounter + threadIdx.x;
            if(my_igp < ngpown)
            {
                int indigp = inv_igp_index[my_igp];
                int igp = indinv[indigp];

                for(int iw = nstart; iw < nend; ++iw)
                {

                    for(int ig = 0; ig < ncouls; ++ig) 
                    { 
                        GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls+igp]);
                        GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                        ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);

                    } //ncouls

                } // iw
            } // ngpown
        }

        for(int iw = nstart; iw < nend; ++iw)
        {
            atomicAdd(&achtemp_re[iw] , achtemp_re_loc[iw] );
            atomicAdd(&achtemp_im[iw] , achtemp_im_loc[iw] );
        }
    }
}

__global__ void NumBandNcouls_kernel( GPUComplex *wtilde_array, GPUComplex *aqsntemp, GPUComplex* aqsmtemp, GPUComplex *I_eps_array, int ncouls, int ngpown, int number_bands, double* wx_array, double* achtemp_re, double* achtemp_im, double* vcoul, int* indinv, int* inv_igp_index, int numThreadsPerBlock, int stride)  
{
    int n1 = blockIdx.x ;

    if(n1 < number_bands)
    {
        int loopOverncouls = 1, leftOverncouls = 0;

        if(ncouls > numThreadsPerBlock)
        {
            loopOverncouls = ncouls / numThreadsPerBlock;
            leftOverncouls = ncouls % numThreadsPerBlock;
        }

        double achtemp_re_loc[nend-nstart];
        double achtemp_im_loc[nend-nstart];
        for(int iw = nstart; iw < nend; ++iw) {achtemp_re_loc[iw] = 0.00; achtemp_im_loc[iw] = 0.00;}

        for(int my_igp = 0; my_igp < ngpown; ++my_igp)
        {
            int indigp = inv_igp_index[my_igp];
            int igp = indinv[indigp];

            for( int x = 0; x < loopOverncouls && threadIdx.x < numThreadsPerBlock; ++x)
            {
                int ig = x*numThreadsPerBlock + threadIdx.x;
                if(ig < ncouls)
                {
                    for(int iw = nstart; iw < nend; ++iw)
                    {
                        GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls+igp]);
                        GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                        ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
                    }
                }
            }
            if(leftOverncouls)
            {
                int ig = loopOverncouls*numThreadsPerBlock + threadIdx.x;
                if(ig < ncouls)
                {
                    for(int iw = nstart; iw < nend; ++iw)
                    {
                        GPUComplex mygpvar1 = d_GPUComplex_conj(aqsmtemp[n1*ncouls+igp]);
                        GPUComplex wdiff = d_doubleMinusGPUComplex(wx_array[iw] , wtilde_array[my_igp*ncouls+ig]);
                        ncoulsKernel(mygpvar1, wdiff, aqsntemp[n1*ncouls+ig], wtilde_array[my_igp*ncouls+ig], I_eps_array[my_igp*ncouls+ig], vcoul[igp], achtemp_re_loc[iw], achtemp_im_loc[iw]);
                    }
                }
            }
        }

        for(int iw = nstart; iw < nend; ++iw)
        {
            atomicAdd(&achtemp_re[iw], achtemp_re_loc[iw]);
            atomicAdd(&achtemp_im[iw], achtemp_im_loc[iw]);
        }
    }
}

void gppKernelGPU( GPUComplex *wtilde_array, GPUComplex *aqsntemp, GPUComplex* aqsmtemp, GPUComplex *I_eps_array, int ncouls, int ngpown, int number_bands, double* wx_array, double *achtemp_re, double *achtemp_im, double *vcoul, int* indinv, int* inv_igp_index, int stride)
{
#if NumBandsKernel
    int numBlocks = number_bands;
    int numThreadsPerBlock = 32;
    printf("Launching a single dimension grid with numBlocks = %d and %d threadsPerBlock \n", numBlocks, numThreadsPerBlock);

    NumberBands_kernel <<< number_bands, numThreadsPerBlock >>> ( wtilde_array, aqsntemp, aqsmtemp, I_eps_array, ncouls, ngpown, number_bands, wx_array, achtemp_re, achtemp_im, vcoul, indinv, inv_igp_index, numThreadsPerBlock, stride);
#endif

#if NgpownKernel
    int numBlocks = ngpown;
    int numThreadsPerBlock = 32;
    printf("Launching a single dimension grid with numBlocks = %d and %d threadsPerBlock \n", numBlocks, numThreadsPerBlock);

    for(int n1 = 0; n1 < number_bands; ++n1)
    {
        cudaNgpown_kernel <<< numBlocks, numThreadsPerBlock >>> ( n1, wtilde_array, aqsntemp, aqsmtemp, I_eps_array, ncouls, ngpown, number_bands, wx_array, achtemp_re, achtemp_im, vcoul, indinv, inv_igp_index, numThreadsPerBlock);
    }
#endif

#if NumBandsNgpownKernel 
    dim3 numBlocks(number_bands, ngpown);
    int numThreadsPerBlock = 32;

    printf("Launching a double dimension grid with numBlocks = (%d, %d) and %d threadsPerBlock \n", number_bands, ngpown, numThreadsPerBlock);
    NumBandNgpown_kernel  <<< numBlocks, numThreadsPerBlock >>> ( wtilde_array, aqsntemp, aqsmtemp, I_eps_array, ncouls, ngpown, number_bands, wx_array, achtemp_re, achtemp_im, vcoul, indinv, inv_igp_index, numThreadsPerBlock, stride);
#endif

#if NgpownNcoulsKernel 
    int numBlocks = number_bands;
//    dim3 numThreadsPerBlock(8, 8);
    int numThreadsPerBlock = 32;

    printf("Launching a single dimension grid with numBlocks =  %d and double dimension threadsPerBlock (8, 8)\n", number_bands);
    NgpownNcouls_kernel  <<< numBlocks, numThreadsPerBlock >>> ( wtilde_array, aqsntemp, aqsmtemp, I_eps_array, ncouls, ngpown, number_bands, wx_array, achtemp_re, achtemp_im, vcoul, indinv, inv_igp_index, numThreadsPerBlock);
#endif

#if NumBandsNcoulsKernel 
    int numBlocks = number_bands;
    int numThreadsPerBlock = 32;

    printf("Launching a double dimension grid with numBlocks = (%d, %d) and %d threadsPerBlock \n", number_bands, ngpown, numThreadsPerBlock);
    NumBandNcouls_kernel  <<< numBlocks, numThreadsPerBlock >>> ( wtilde_array, aqsntemp, aqsmtemp, I_eps_array, ncouls, ngpown, number_bands, wx_array, achtemp_re, achtemp_im, vcoul, indinv, inv_igp_index, numThreadsPerBlock, stride);
#endif
}

void till_nvbandKernel(GPUComplex *aqsmtemp, GPUComplex *aqsntemp, GPUComplex *asxtemp, int *inv_igp_index, int *indinv, GPUComplex *wtilde_array, double *wx_array, GPUComplex *I_eps_array, int ncouls, int nvband, int ngpown, double *vcoul)
{
    dim3 numBlocks(nvband, ngpown);
    int numThreadsPerBlock = ncouls;

    d_flagOCC_solver<<< numBlocks, numThreadsPerBlock>>>(wx_array, wtilde_array, asxtemp, aqsmtemp, aqsntemp, I_eps_array, inv_igp_index, indinv, ncouls, nvband, ngpown, vcoul);
}
