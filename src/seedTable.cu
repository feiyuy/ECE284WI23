#include "seedTable.cuh"
#include <stdio.h>
#include <math.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>


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
    int* d_seq1,
    int* d_seq2,
    int *mat,
    int L1,
    int L2
    //L1<=L2
) {

    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;

    int x, y, up, left, diagonal, temp, result;

    for (int i=1; i<L2; i++){
        if (thread_index < (i < L1-1? i:L1-1)){
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

            temp = (up < left? up:left);
            result = (temp < diagonal? temp:diagonal);
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

            temp = (up < left? up:left);
            result = (temp < diagonal? temp:diagonal);
            mat[y*L1+x] = result;
        }

        __syncthreads();
    }
}


void compute(
    int *seq1,
    int *seq2,
    int l1,
    int l2,
    int *mat
) {

    int *d_mat, *d_seq1, *d_seq2;
    int L1, L2;
    cudaMalloc(&d_mat, l1*l2*sizeof(int));
    if (l1<l2){
        L1 = l1;
        L2 = l2;
        cudaMalloc(&d_seq1, L1*sizeof(int));
        cudaMalloc(&d_seq2, L2*sizeof(int));
        cudaMemcpy(d_seq1, seq1, L1*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_seq2, seq2, L2*sizeof(int), cudaMemcpyHostToDevice);
    }
    else{
        L1 = l2;
        L2 = l1;
        cudaMalloc(&d_seq1, L1*sizeof(int));
        cudaMalloc(&d_seq2, L2*sizeof(int));
        cudaMemcpy(d_seq1, seq2, L1*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_seq2, seq1, L2*sizeof(int), cudaMemcpyHostToDevice);
    }

    initRows<<<5, L1/5>>>(d_mat, L1);
    initCols<<<5, L2/5>>>(d_mat, L1, L2);
    matrixFill<<<5, L2/5>>>(d_seq1, d_seq2, d_mat, L1, L2);

    cudaMemcpy(mat, d_mat, l1*l2*sizeof(int), cudaMemcpyDeviceToHost);
}