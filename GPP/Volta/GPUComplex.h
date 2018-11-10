#ifndef __GPUCOMPLEX
#define __GPUCOMPLEX

#include <iostream>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <cmath>
#include <complex>
#include <omp.h>
#include <ctime>
#include <chrono>


#include <vector_types.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

#define NumBandsKernel 0
#define NgpownKernel 0
#define NumBandsNgpownKernel 1
#define NgpownNcoulsKernel 0
#define NumBandsNcoulsKernel 0

#define nstart 0
#define nend 6

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
        file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif

    return;
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
        file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }

    // More careful checking. However, this will affect performance.
    // Comment away if needed. - Rahul - commented the below deviceSynchronize
//    err = cudaDeviceSynchronize();
    if( cudaSuccess != err )
    {
        fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif
    return;
}


class GPUComplex : public double2{

    private : 

public:
//    double2 complNum;

explicit GPUComplex () {
    x = 0.00;
    y = 0.00;
}


__host__ __device__ explicit GPUComplex(const double& a, const double& b) {
    x = a;
    y = b;
}

GPUComplex(const GPUComplex& src) {
    x = src.x;
    y = src.y;
}

GPUComplex& operator =(const GPUComplex& src) {
    x = src.x;
    y = src.y;

    return *this;
}

GPUComplex& operator +=(const GPUComplex& src) {
    x = src.x + this->x;
    y = src.y + this->y;

    return *this;
}

GPUComplex& operator -=(const GPUComplex& src) {
    x = src.x - this->x;
    y = src.y - this->y;

    return *this;
}

GPUComplex& operator -() {
    x = -this->x;
    y = -this->y;

    return *this;
}

GPUComplex& operator ~() {
    return *this;
}

void print() const {
    printf("( %f, %f) ", this->x, this->y);
    printf("\n");
}

double abs(const GPUComplex& src) {

    double re_this = src.x * src.x;
    double im_this = src.y * src.y;

    double result = (re_this+im_this);
    return result;
}

double get_real() const
{
    return this->x;
}

double get_imag() const
{
    return this->y;
}

void set_real(double val)
{
    this->x = val;
}

void set_imag(double val) 
{
    this->y = val;
}

    friend inline GPUComplex GPUComplex_square(GPUComplex& src) ;
    friend inline GPUComplex GPUComplex_conj(const GPUComplex& src) ;
    friend inline GPUComplex GPUComplex_product(const GPUComplex& a, const GPUComplex& b) ;
    friend inline double GPUComplex_abs(const GPUComplex& src) ;
    friend inline GPUComplex GPUComplex_mult(GPUComplex& a, double b, double c) ;
    friend inline GPUComplex GPUComplex_mult(const GPUComplex& a, double b) ;
    friend inline void GPUComplex_fma(GPUComplex& a, const GPUComplex& b, const GPUComplex& c) ;
    friend inline void GPUComplex_fms(GPUComplex& a, const GPUComplex& b, const GPUComplex& c) ;
    friend inline GPUComplex doubleMinusGPUComplex(const double &a, GPUComplex& src) ;
    friend inline GPUComplex doublePlusGPUComplex(double a, GPUComplex& src) ;
    friend inline double GPUComplex_real( const GPUComplex& src) ;
    friend inline double GPUComplex_imag( const GPUComplex& src) ;

    
//Device Functions 
    friend __device__ GPUComplex d_GPUComplex_square(GPUComplex& src) ;
    friend __device__ GPUComplex d_GPUComplex_conj(const GPUComplex& src) ;
    friend __device__ GPUComplex d_GPUComplex_product(const GPUComplex& a, const GPUComplex& b) ;
    friend __device__ double d_GPUComplex_abs(const GPUComplex& src) ;
    friend __device__ GPUComplex d_GPUComplex_mult(GPUComplex& a, double b, double c) ;
    friend __device__ GPUComplex d_GPUComplex_mult(const GPUComplex& a, double b) ;
    friend __device__ void d_GPUComplex_fma(GPUComplex& a, const GPUComplex& b, const GPUComplex& c) ;
    friend __device__ void d_GPUComplex_fms(GPUComplex& a, const GPUComplex& b, const GPUComplex& c) ;
    friend __device__ GPUComplex d_doubleMinusGPUComplex(const double &a, GPUComplex& src) ;
    friend __device__ GPUComplex d_doublePlusGPUComplex(double a, GPUComplex& src) ;
    friend __device__ double d_GPUComplex_real( const GPUComplex& src) ;
    friend __device__ double d_GPUComplex_imag( const GPUComplex& src) ;
    friend __device__ void d_GPUComplex_plusEquals( GPUComplex& a, const GPUComplex & b); 
    friend __device__ void d_GPUComplex_Equals( GPUComplex& a, const GPUComplex & b); 
    friend __device__ void d_print( const GPUComplex& src) ;
    friend __device__ void ncoulsKernel(GPUComplex& mygpvar1, GPUComplex& wdiff, GPUComplex& aqsntemp_index, GPUComplex& wtilde_array_index, GPUComplex& I_eps_array_index, double vcoul_igp, double& achtemp_re_loc, double& achtemp_im_loc);
    friend __device__ void ncoulsKernel(GPUComplex& mygpvar1, double vcoul_igp, double& achtemp_re_loc, double& achtemp_im_loc);

};
//Inline functions have to be defined in the same file as the declaration

/*
 * Return the square of a complex number 
 */
GPUComplex GPUComplex_square(GPUComplex& src) {
    double re_this = src.x ;
    double im_this = src.y ;

    GPUComplex result(re_this*re_this - im_this*im_this, 2*re_this*im_this);

    return result;
}

/*
 * Return the conjugate of a complex number 
 */
GPUComplex GPUComplex_conj(const GPUComplex& src) {

    double re_this = src.x;
    double im_this = -1 * src.y;

    GPUComplex result(re_this, im_this);
    return result;
}


/*
 * Return the product of 2 complex numbers 
 */
GPUComplex GPUComplex_product(const GPUComplex& a, const GPUComplex& b) {

    double re_this = a.x * b.x - a.y*b.y ;
    double im_this = a.x * b.y + a.y*b.x ;

    GPUComplex result(re_this, im_this);
    return result;
}

/*
 * Return the absolute of a complex number 
 */
double GPUComplex_abs(const GPUComplex& src) {
    double re_this = src.x * src.x;
    double im_this = src.y * src.y;

    double result = (re_this+im_this);
    return result;
}

/*
 *  result = a * b * c (a = complex ; b,c = double) 
 */
GPUComplex GPUComplex_mult(GPUComplex& a, double b, double c) {

    GPUComplex result(a.x * b * c, a.y * b * c);
    return result;

}

/*
 * Return the complex number c = a * b (a is complex, b is double) 
 */
GPUComplex GPUComplex_mult(const GPUComplex& a, double b) {

   GPUComplex result(a.x*b, a.y*b);
   return result;

}

/*
 * Return the complex number a += b * c  
 */
void GPUComplex_fma(GPUComplex& a, const GPUComplex& b, const GPUComplex& c) {
    double re_this = b.x * c.x - b.y*c.y ;
    double im_this = b.x * c.y + b.y*c.x ;

    GPUComplex mult_result(re_this, im_this);

    a.x += mult_result.x;
    a.y += mult_result.y;
}

/*
 * Return the complex number a -= b * c  
 */
void GPUComplex_fms(GPUComplex& a, const GPUComplex& b, const GPUComplex& c) {
    double re_this = b.x * c.x - b.y*c.y ;
    double im_this = b.x * c.y + b.y*c.x ;

    GPUComplex mult_result(re_this, im_this);

    a.x -= mult_result.x;
    a.y -= mult_result.y;
}


GPUComplex doubleMinusGPUComplex(const double &a, GPUComplex& src) {
    GPUComplex result(a - src.x, 0 - src.y);
    return result;
}

GPUComplex doublePlusGPUComplex(double a, GPUComplex& src) {
    GPUComplex result(a + src.x, 0 + src.y);
    return result;
}

double GPUComplex_real( const GPUComplex& src) {
    return src.x;
}

double GPUComplex_imag( const GPUComplex& src) {
    return src.y;
}

void gppKernelGPU( GPUComplex *wtilde_array, GPUComplex *aqsntemp, GPUComplex* aqsmtemp, GPUComplex *I_eps_array, int ncouls, int ngpown, int number_bands, double* wx_array, double *achtemp_re, double *achtemp_im, double *vcoul, int* indinv, int* inv_igp_index, int stride);

void till_nvbandKernel(GPUComplex *aqsmtemp, GPUComplex *aqsntemp, GPUComplex *asxtemp, int *inv_igp_index, int *indinv, GPUComplex *wtilde_array, double *wx_array, GPUComplex *I_eps_array, int ncouls, int nvband, int ngpown, double *d_vcoul);

#endif
