#ifndef KERNEL1_H
#define KERNEL1_H

#ifdef ERT_GPU
extern int gpu_blocks;
extern int gpu_threads;
#endif

#define KERNEL1(a,b,c)   ((a) = (b) + (c))
#define KERNEL2(a,b,c)   ((a) = (a)*(b) + (c))

void initialize(uint64_t nsize,
                double* __restrict__ array,
                double value);

#ifdef ERT_GPU
void gpuKernel(uint64_t nsize,
               uint64_t ntrials,
               double* __restrict__ array,
               int* bytes_per_elem,
               int* mem_accesses_per_elem);
#else
void kernel(uint64_t nsize,
            uint64_t ntrials,
            double* __restrict__ array,
            int* bytes_per_elem,
            int* mem_accesses_per_elem);
#endif

#endif
