#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#define MIN(x,y) ((x) < (y) ? (x) : (y))

__global__ void initRows(
    int *mat,
    int L1
) {

    int total_threads = gridDim.x * blockDim.x;

    int entry_per_thread = L1 / total_threads + 1;

    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;

    for (int i=entry_per_thread*thread_index; i<entry_per_thread*(thread_index+1); i++){
        if (i<L1){
            mat[i] = i;
        }
    }
}

__global__ void initCols(
    int *mat,
    int L1,
    int L2
){
    
    int total_threads = gridDim.x * blockDim.x;

    int entry_per_thread = L2 / total_threads + 1;

    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;

    for (int i=entry_per_thread*thread_index; i<entry_per_thread*(thread_index+1); i++){
        if (i<L2){
            mat[i*L1] = i;
        }
    }
}

__global__ void matrixFill(
    unit32_t* d_seq1,
    unit32_t* d_seq2,
    int *mat,
    int L1,
    int L2
    //L1<=L2
) {

    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;

    for (int i=1; i<L2; i++){
        if (thread_index < MIN(i, L1-1)){
            x = 1 + thread_index;
            y = i - thread_index;

            up = mat[(y-1)*L1+x];
            left = mat[y*L1+x-1];
            if (d_seq1[x] == d_seq2[y]){
                diagonal = mat[(y-1)*L1+x-1];
            }
            else{
                diagonal = mat[(y-1)*L1+x-1] + 1;
            }

            temp = MIN(up, left);
            result = MIN(temp, diagonal);
            mat[y*L1+x] = result;
        }

        __syncthreads();
    }

    for (int j=1; j<L1-1; j++){
        if (thread_index < j){
            x = 1 + L1 - j + thread_index;
            y = L2 - 1 - thread_index;

            up = mat[(y-1)*L1+x];
            left = mat[y*L1+x-1];
            if (d_seq1[x] == d_seq2[y]){
                diagonal = mat[(y-1)*L1+x-1];
            }
            else{
                diagonal = mat[(y-1)*L1+x-1] + 1;
            }

            temp = MIN(up, left);
            result = MIN(temp, diagonal);
            mat[y*L1+x] = result;
        }

        __syncthreads();
    }
}
