#include <iostream>
#include <cstdlib>
#include <memory>

#include <iomanip>
#include <cmath>
#include <complex>
#include <omp.h>
#include <ctime>
#include <chrono>
#include<inttypes.h>

#include "GPUComplex.h"

using namespace std;

GPUComplex** allocateMemGPU(size_t size)
{
    void **d_src;
    if(cudaMalloc((void**) &d_src, size) != cudaSuccess)
    {
        return NULL;
    }

    return (GPUComplex**) d_src;
}

inline void reduce_achstemp(int n1, int number_bands, int* inv_igp_index, int ncouls, GPUComplex  *aqsmtemp, GPUComplex *aqsntemp, GPUComplex *I_eps_array, GPUComplex achstemp,  int* indinv, int ngpown, double* vcoul)
{
    double to1 = 1e-6;
    GPUComplex schstemp(0.0, 0.0);;

    for(int my_igp = 0; my_igp< ngpown; my_igp++)
    {
        GPUComplex schs(0.0, 0.0);
        GPUComplex matngmatmgp(0.0, 0.0);
        GPUComplex matngpmatmg(0.0, 0.0);
        GPUComplex halfinvwtilde, delw, ssx, sch, wdiff, cden , eden, mygpvar1, mygpvar2;
        int indigp = inv_igp_index[my_igp];
        int igp = indinv[indigp];
        if(indigp == ncouls)
            igp = ncouls-1;

        if(!(igp > ncouls || igp < 0)){

            GPUComplex mygpvar2, mygpvar1;
            mygpvar1 = GPUComplex_conj(aqsmtemp[n1*ncouls+igp]);
            mygpvar2 = aqsntemp[n1*ncouls+igp];



            schs = I_eps_array[my_igp*ncouls+igp];
            matngmatmgp = GPUComplex_product(mygpvar1, aqsntemp[n1*ncouls+igp]);


            if(GPUComplex_abs(schs) > to1)
                GPUComplex_fma(schstemp, matngmatmgp, schs);
            }
            else 
            {
                for(int ig=1; ig<ncouls; ++ig)
                {
                    GPUComplex mult_result(GPUComplex_product(I_eps_array[my_igp*ncouls+ig] , mygpvar1));
                    GPUComplex_fms(schstemp,aqsntemp[n1*ncouls+igp], mult_result); 
                }
            }

        schstemp = GPUComplex_mult(schstemp, vcoul[igp], 0.5);
        achstemp += schstemp;
    }
}

inline void flagOCC_solver(double wxt, GPUComplex *wtilde_array, int my_igp, int n1, GPUComplex *aqsmtemp, GPUComplex *aqsntemp, GPUComplex *I_eps_array, GPUComplex &ssxt, GPUComplex &scht,int ncouls, int igp, int number_bands, int ngpown)
{
    GPUComplex expr0(0.00, 0.00);
    GPUComplex expr(0.5, 0.5);
    GPUComplex matngmatmgp(0.0, 0.0);
    GPUComplex matngpmatmg(0.0, 0.0);

    for(int ig=0; ig<ncouls; ++ig)
    {
        GPUComplex wtilde = wtilde_array[my_igp*ncouls+ig];
        GPUComplex wtilde2 = GPUComplex_square(wtilde);
        GPUComplex Omega2 = GPUComplex_product(wtilde2,I_eps_array[my_igp*ncouls+ig]);
        GPUComplex mygpvar1 = GPUComplex_conj(aqsmtemp[n1*ncouls+igp]);
        GPUComplex mygpvar2 = aqsmtemp[n1*ncouls+igp];
        GPUComplex matngmatmgp = GPUComplex_product(aqsntemp[n1*ncouls+ig] , mygpvar1);
        if(ig != igp) matngpmatmg = GPUComplex_product(GPUComplex_conj(aqsmtemp[n1*ncouls+ig]) , mygpvar2);

        double delw2, scha_mult, ssxcutoff;
        double to1 = 1e-6;
        double sexcut = 4.0;
        double gamma = 0.5;
        double limitone = 1.0/(to1*4.0);
        double limittwo = pow(0.5,2);
        GPUComplex sch, ssx;
    
        GPUComplex wdiff = doubleMinusGPUComplex(wxt , wtilde);
    
        GPUComplex cden = wdiff;
        double rden = 1/GPUComplex_real(GPUComplex_product(cden , GPUComplex_conj(cden)));
        GPUComplex delw = GPUComplex_mult(GPUComplex_product(wtilde , GPUComplex_conj(cden)) , rden);
        double delwr = GPUComplex_real(GPUComplex_product(delw , GPUComplex_conj(delw)));
        double wdiffr = GPUComplex_real(GPUComplex_product(wdiff , GPUComplex_conj(wdiff)));
    
        if((wdiffr > limittwo) && (delwr < limitone))
        {
            sch = GPUComplex_product(delw , I_eps_array[my_igp*ngpown+ig]);
            double cden = std::pow(wxt,2);
            rden = std::pow(cden,2);
            rden = 1.00 / rden;
            ssx = GPUComplex_mult(Omega2 , cden , rden);
        }
        else if (delwr > to1)
        {
            sch = expr0;
            cden = GPUComplex_mult(GPUComplex_product(wtilde2, doublePlusGPUComplex((double)0.50, delw)), 4.00);
            rden = GPUComplex_real(GPUComplex_product(cden , GPUComplex_conj(cden)));
            rden = 1.00/rden;
            ssx = GPUComplex_product(GPUComplex_product(-Omega2 , GPUComplex_conj(cden)), GPUComplex_mult(delw, rden));
        }
        else
        {
            sch = expr0;
            ssx = expr0;
        }
    
        ssxcutoff = GPUComplex_abs(I_eps_array[my_igp*ngpown+ig]) * sexcut;
        if((GPUComplex_abs(ssx) > ssxcutoff) && (wxt < 0.00)) ssx = expr0;

        ssxt += GPUComplex_product(matngmatmgp , ssx);
        scht += GPUComplex_product(matngmatmgp , sch);
    }
}

void gppKernelCPU( GPUComplex *wtilde_array, GPUComplex *aqsntemp, GPUComplex *I_eps_array, int ncouls, double wxt, double& achtemp_re_iw, double& achtemp_im_iw, int my_igp, GPUComplex mygpvar1, int n1, double vcoul_igp)
{
    GPUComplex scht(0.00, 0.00);
    for(int ig = 0; ig<ncouls; ++ig)
    {

        GPUComplex wdiff = doubleMinusGPUComplex(wxt , wtilde_array[my_igp*ncouls+ig]);
        double rden = GPUComplex_real(GPUComplex_product(wdiff, GPUComplex_conj(wdiff)));
        rden = 1/rden;
        GPUComplex delw = GPUComplex_mult(GPUComplex_product(wtilde_array[my_igp*ncouls+ig] , GPUComplex_conj(wdiff)), rden); 
        
        scht += GPUComplex_mult(GPUComplex_product(GPUComplex_product(mygpvar1 , aqsntemp[n1*ncouls+ig]), GPUComplex_product(delw , I_eps_array[my_igp*ncouls+ig])), 0.5);
    }
    achtemp_re_iw += GPUComplex_real( GPUComplex_mult(scht , vcoul_igp));
    achtemp_im_iw += GPUComplex_imag( GPUComplex_mult(scht , vcoul_igp));
}

int main(int argc, char** argv)
{

    if (argc != 6)
    {
        std::cout << "The correct form of input is : " << endl;
        std::cout << " ./a.out <number_bands> <number_valence_bands> <number_plane_waves> <nodes_per_mpi_group> <stride> " << endl;
        exit (0);
    }

    printf("********Executing Cuda version of the Kernel*********\n");

    auto start_totalTime = std::chrono::high_resolution_clock::now();
    const int number_bands = atoi(argv[1]);
    const int nvband = atoi(argv[2]);
    const int ncouls = atoi(argv[3]);
    const int nodes_per_group = atoi(argv[4]);
    const int stride = atoi(argv[5]);
    const int npes = 1; 
    const int ngpown = ncouls / (nodes_per_group * npes); 
    const double e_lk = 10;
    const double dw = 1;

    double to1 = 1e-6, \
    gamma = 0.5, \
    sexcut = 4.0;
    double limitone = 1.0/(to1*4.0), \
    limittwo = pow(0.5,2);
    const double e_n1kq= 6.0; 

    //Printing out the params passed.
    std::cout << "number_bands = " << number_bands \
        << "\t nvband = " << nvband \
        << "\t ncouls = " << ncouls \
        << "\t nodes_per_group  = " << nodes_per_group \
        << "\t ngpown = " << ngpown \
        << "\t nend = " << nend \
        << "\t nstart = " << nstart \
        << "\t gamma = " << gamma \
        << "\t sexcut = " << sexcut \
        << "\t limitone = " << limitone \
        << "\t limittwo = " << limittwo << endl;


    //ALLOCATE statements from fortran gppkernel.
    
   
    GPUComplex expr0(0.00, 0.00);
    GPUComplex expr(0.5, 0.5);

    int *inv_igp_index = new int[ngpown];
    int *indinv = new int[ncouls+1];

    GPUComplex *acht_n1_loc = new GPUComplex[number_bands];
    GPUComplex *achtemp = new GPUComplex[(nend-nstart)];
    GPUComplex *aqsmtemp = new GPUComplex[number_bands*ncouls];
    GPUComplex *aqsntemp = new GPUComplex[number_bands*ncouls];
    GPUComplex *I_eps_array = new GPUComplex[ngpown*ncouls];
    GPUComplex *wtilde_array = new GPUComplex[ngpown*ncouls];
    GPUComplex *ssx_array = new GPUComplex[(nend-nstart)];
    GPUComplex *ssxa = new GPUComplex[ncouls];
    GPUComplex achstemp;

    double *achtemp_re = new double[(nend-nstart)];
    double *wx_array = new double[(nend-nstart)];
    double *achtemp_im = new double[(nend-nstart)];
    double *vcoul = new double[ncouls];

    printf("Executing CUDA version of the Kernel stride = %d\n", stride);
//Data Structures on Device
    GPUComplex *d_wtilde_array, *d_aqsntemp, *d_aqsmtemp, *d_I_eps_array, *d_asxtemp;
    double *d_achtemp_re, *d_achtemp_im, *d_vcoul, *d_wx_array;
    int *d_inv_igp_index, *d_indinv;

    CudaSafeCall(cudaMalloc((void**) &d_wtilde_array, ngpown*ncouls*sizeof(GPUComplex)));
    CudaSafeCall(cudaMalloc((void**) &d_I_eps_array, ngpown*ncouls*sizeof(GPUComplex)));
    CudaSafeCall(cudaMalloc((void**) &d_aqsntemp, number_bands*ncouls*sizeof(GPUComplex)));
    CudaSafeCall(cudaMalloc((void**) &d_aqsmtemp, number_bands*ncouls*sizeof(GPUComplex)));
    CudaSafeCall(cudaMalloc((void**) &d_achtemp_re, (nend-nstart)*sizeof(double)));
    CudaSafeCall(cudaMalloc((void**) &d_achtemp_im, (nend-nstart)*sizeof(double)));
    CudaSafeCall(cudaMalloc((void**) &d_wx_array, (nend-nstart)*sizeof(double)));
    CudaSafeCall(cudaMalloc((void**) &d_vcoul, ncouls*sizeof(double)));
    CudaSafeCall(cudaMalloc((void**) &d_indinv, (ncouls+1)*sizeof(int)));
    CudaSafeCall(cudaMalloc((void**) &d_inv_igp_index, ngpown*sizeof(int)));
    CudaSafeCall(cudaMalloc((void**) &d_asxtemp, (nend-nstart)*sizeof(double)));

    double occ=1.0;
    bool flag_occ;
    double achstemp_real = 0.00, achstemp_imag = 0.00;
//    cout << "Size of wtilde_array = " << (ncouls*ngpown*2.0*8) / pow(1024,2) << " Mbytes" << endl;
//    cout << "Size of aqsntemp = " << (ncouls*number_bands*2.0*8) / pow(1024,2) << " Mbytes" << endl;
//    cout << "Size of I_eps_array array = " << (ncouls*ngpown*2.0*8) / pow(1024,2) << " Mbytes" << endl;
//
   for(int i=0; i<number_bands; i++)
       for(int j=0; j<ncouls; j++)
       {
           aqsmtemp[i*ncouls+j].x = 0.5;
           aqsmtemp[i*ncouls+j].y = 0.5;
           aqsntemp[i*ncouls+j].x = 0.5;
           aqsntemp[i*ncouls+j].y = 0.5;
       }

   for(int i=0; i<ngpown; i++)
       for(int j=0; j<ncouls; j++)
       {
           I_eps_array[i*ncouls+j].x = 0.5;
           I_eps_array[i*ncouls+j].y = 0.5;
           wtilde_array[i*ncouls+j].x = 0.5;
           wtilde_array[i*ncouls+j].y = 0.5;
       }

   for(int i=0; i<ncouls; i++)
       vcoul[i] = 1.0;


    for(int ig=0, tmp=1; ig < ngpown; ++ig,tmp++)
        inv_igp_index[ig] = (ig+1) * ncouls / ngpown;

    //Do not know yet what this array represents
    for(int ig=0 ; ig<ncouls; ++ig)
        indinv[ig] = ig;
        indinv[ncouls] = ncouls-1;

    for(int iw=nstart; iw<nend; ++iw)
    {
        wx_array[iw] = e_lk - e_n1kq + dw*((iw+1)-2);
        if(wx_array[iw] < to1) wx_array[iw] = to1;
    }

    auto start_withDataMovement = std::chrono::high_resolution_clock::now();
    float mem_alloc = 0.00;

//Start memcpyToDevice 

    CudaSafeCall(cudaMemcpy(d_wtilde_array, wtilde_array, ngpown*ncouls*sizeof(GPUComplex), cudaMemcpyHostToDevice));

    CudaSafeCall(cudaMemcpy(d_I_eps_array, I_eps_array, ngpown*ncouls*sizeof(GPUComplex), cudaMemcpyHostToDevice));
    mem_alloc += 2*ngpown*ncouls*sizeof(GPUComplex);

    CudaSafeCall(cudaMemcpy(d_aqsmtemp, aqsmtemp, number_bands*ncouls*sizeof(GPUComplex), cudaMemcpyHostToDevice));

    CudaSafeCall(cudaMemcpy(d_aqsntemp, aqsntemp, number_bands*ncouls*sizeof(GPUComplex), cudaMemcpyHostToDevice));
    mem_alloc += 2*number_bands*ncouls*sizeof(GPUComplex);

    CudaSafeCall(cudaMemcpy(d_indinv, indinv, (ncouls+1)*sizeof(int), cudaMemcpyHostToDevice));
    mem_alloc += ncouls*sizeof(int);

    CudaSafeCall(cudaMemcpy(d_inv_igp_index, inv_igp_index, ngpown*sizeof(int), cudaMemcpyHostToDevice));
    mem_alloc += ngpown*sizeof(int);

    CudaSafeCall(cudaMemcpy(d_vcoul, vcoul, ncouls*sizeof(double), cudaMemcpyHostToDevice));
    mem_alloc += ncouls*sizeof(double);

    CudaSafeCall(cudaMemcpy(d_wx_array, wx_array, (nend-nstart)*sizeof(double), cudaMemcpyHostToDevice));

    CudaSafeCall(cudaMemcpy(d_achtemp_re, achtemp_re, (nend-nstart)*sizeof(double), cudaMemcpyHostToDevice));

    CudaSafeCall(cudaMemcpy(d_achtemp_im, achtemp_im, (nend-nstart)*sizeof(double), cudaMemcpyHostToDevice));
    mem_alloc += 3*3*sizeof(double);

    mem_alloc /= (1024*1024*1024);

    printf("mem_alloc = %f GBs\n", mem_alloc);
//Start Kernel 
    auto start_kernelTiming = std::chrono::high_resolution_clock::now();

    till_nvbandKernel(d_aqsmtemp, d_aqsntemp, d_asxtemp, d_inv_igp_index, d_indinv, d_wtilde_array, d_wx_array, d_I_eps_array, ncouls, nvband, ngpown, d_vcoul);

    gppKernelGPU( d_wtilde_array, d_aqsntemp, d_aqsmtemp, d_I_eps_array, ncouls, ngpown, number_bands, d_wx_array, d_achtemp_re, d_achtemp_im, d_vcoul, d_indinv, d_inv_igp_index, stride);

    cudaDeviceSynchronize();
    std::chrono::duration<double> elapsed_kernelTiming = std::chrono::high_resolution_clock::now() - start_kernelTiming;

//Start memcpyToHost 
    CudaSafeCall(cudaMemcpy(achtemp_im, d_achtemp_im, (nend-nstart)*sizeof(double), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(achtemp_re, d_achtemp_re, (nend-nstart)*sizeof(double), cudaMemcpyDeviceToHost));

    printf(" \n Cuda Kernel Final achtemp\n");
    for(int iw=nstart; iw<nend; ++iw)
    {
        achtemp[iw] = GPUComplex(achtemp_re[iw], achtemp_im[iw]);
//        achtemp[iw].print();
    }
    achtemp[0].print();

    std::chrono::duration<double> elapsed_totalTime = std::chrono::high_resolution_clock::now() - start_totalTime;

    cout << "********** Kernel Time Taken **********= " << elapsed_kernelTiming.count() << " secs" << endl;
    cout << "********** Total Time Taken **********= " << elapsed_totalTime.count() << " secs" << endl;

    cudaFree(d_wtilde_array);
    cudaFree(d_aqsntemp);
    cudaFree(d_aqsntemp);
    cudaFree(d_asxtemp);
    cudaFree(d_I_eps_array);
    cudaFree(d_achtemp_re);
    cudaFree(d_achtemp_im);
    cudaFree(d_vcoul);
    cudaFree(d_inv_igp_index);
    cudaFree(d_indinv);

    free(acht_n1_loc);
    free(achtemp);
    free(aqsmtemp);
    free(aqsntemp);
    free(I_eps_array);
    free(wtilde_array);
    free(vcoul);
    free(ssx_array);

    return 0;
}
