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
    Timer timer;
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
    
    timer.Start();
    fprintf(stdout, "\nCompute Levinstein distance in GPU.\n");
    initRows<<<5, L1/5>>>(d_mat, L1);
    initCols<<<5, L2/5>>>(d_mat, L1, L2);
    matrixFill<<<5, L2/5>>>(d_seq1, d_seq2, d_mat, L1, L2);
    fprintf(stdout, "\nCompleted in %ld usec \n\n", timer.Stop());

    cudaMemcpy(mat, d_mat, l1*l2*sizeof(int), cudaMemcpyDeviceToHost);
}


__global__ void initBoundary(
    int *mat
) {

    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    mat[thread_index] = thread_index;
}

__global__ void boundaryFill(
    int* d_seq1,
    int* d_seq2,
    int *d_edge_up1,
    int *d_edge_up2,
    int *d_edge_left1,
    int *d_edge_left2,
    int *d_buffer1,
    int *d_buffer2,
    int *d_buffer3
) {

    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;

    int x, y, up, left, diagonal, temp, result;

    for (int i=1; i<8; i++){
        if (thread_index < i){

            d_buffer3 = d_buffer2;
            d_buffer2 = d_buffer1;

            x = 1 + thread_index;
            y = i - thread_index;

            if (x == 1 && y == 1){
                up = d_edge_up1[1];
                left = d_edge_left1[1];
                diagonal = d_edge_up1[0];
            }
            else if (x == 1){
                up = d_buffer2[0];
                left = d_edge_left1[i];
                diagonal = d_edge_left1[i-1];
            }
            else if (y == 1){
                up = d_edge_up1[i];
                left = d_buffer2[i-2];
                diagonal = d_edge_up1[i-1];
            }
            else{
                up = d_buffer2[thread_index];
                left = d_buffer2[thread_index-1];
                diagonal = d_buffer3[thread_index-1];
            }
            if (d_seq1[x] != d_seq2[y]){
                diagonal = diagonal + 1;
            }

            temp = (up < left? up:left);
            result = (temp < diagonal? temp:diagonal);
            d_buffer1[thread_index] = result;

            if (x == 4){
                d_edge_left2[y] = result;
            }
            if (y == 4){
                d_edge_up2[x] == result;
            }
        }

        __syncthreads();
    }

    for (int j=1; j<7; j++){
        if (thread_index < j){

            d_buffer3 = d_buffer2;
            d_buffer2 = d_buffer1;

            x = 1 + j + thread_index;
            y = 7 - thread_index;

            up = d_buffer2[thread_index+1];
            left = d_buffer2[thread_index];
            diagonal = d_buffer3[thread_index+1];

            if (d_seq1[x] != d_seq2[y]){
                diagonal = diagonal + 1;
            }

            temp = (up < left? up:left);
            result = (temp < diagonal? temp:diagonal);
            d_buffer1[thread_index] = result;

            if (x == 4){
                d_edge_left2[y] = result;
            }
            if (y == 4){
                d_edge_up2[x] == result;
            }
        }

        __syncthreads();
    }
}

__global__ void initRow(
    int *d_edge_up,
    int *d_mat
) {
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    d_mat[thread_index] = d_edge_up[thread_index];
}

__global__ void initCol(
    int *d_edge_left,
    int *d_mat
){

    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;

    d_mat[thread_index*4] = d_edge_left[thread_index];
}

void seg_compute(
    int *seq1,
    int *seq2,
    char *track
) {
    int *d_edge_up1, *d_edge_up2, *d_edge_left1, *d_edge_left2;
    cudaMalloc(&d_edge_up1, 8*sizeof(int));
    cudaMalloc(&d_edge_up2, 8*sizeof(int));
    cudaMalloc(&d_edge_left1, 8*sizeof(int));
    cudaMalloc(&d_edge_left2, 8*sizeof(int));

    int *d_buffer1, *d_buffer2, *d_buffer3;
    cudaMalloc(&d_buffer1, 8*sizeof(int));
    cudaMalloc(&d_buffer2, 8*sizeof(int));
    cudaMalloc(&d_buffer3, 8*sizeof(int));

    int *d_seq1, *d_seq2;
    cudaMalloc(&d_seq1, 8*sizeof(int));
    cudaMalloc(&d_seq2, 8*sizeof(int));
    cudaMemcpy(d_seq1, seq1, 8*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, seq2, 8*sizeof(int), cudaMemcpyHostToDevice);

    initBoundary<<<2, 4>>>(d_edge_up1);
    initBoundary<<<2, 4>>>(d_edge_left1);
    boundaryFill<<<2, 4>>>(d_seq1, d_seq2, d_edge_up1, d_edge_up2, d_edge_left1, d_edge_left2, d_buffer1, d_buffer2, d_buffer3);

    int *d_mat, *mat;
    cudaMalloc(&d_mat, 4*4*sizeof(int));

    initCol<<<2, 2>>>(d_edge_left2+4, d_mat);
    initRow<<<2, 2>>>(d_edge_up2+4, d_mat);

    matrixFill<<<1, 3>>>(d_seq1+4, d_seq2+4, d_mat, 4, 4);

    cudaMemcpy(mat, d_mat, 4*4*sizeof(int), cudaMemcpyDeviceToHost);

    //initial cordinate
    int x = 3;
    int y = 3;
    int temp_x, temp_y, temp_value;
    char temp_direction;

    int index = 0;

    while (x!=0 and y!=0){
        if (mat[(y-1)*4+x] > mat[y*4+x-1]){
            temp_x = x-1;
            temp_y = y;
            temp_value = mat[y*4+x-1];
            temp_direction = 'U'; //up
        }
        else{
            temp_x = x;
            temp_y = y-1;
            temp_value = mat[(y-1)*4+x];
            temp_direction = 'L'; //left
        }
        if (mat[(y-1)*4+x-1] > temp_value){
            x = temp_x;
            y = temp_y;
            track[index] = temp_direction;
        }
        else{
            x = x-1;
            y = y-1;
            track[index] = 'D';
        }
        index = index+1;
    }

    int next;

    if (x == 0 and y ==0){
        if (track[index-1] == 'D'){
            next = 0;
        }
        else if (track[index-1] == 'U'){
            next = 1;
        }
        else{
            next = 2;
        }

    }
    else if (x == 0){
        next = 1;
    }
    else{
        next = 2;
    }

    if (next == 0){
        initCol<<<2, 2>>>(d_edge_left2, d_mat);
        initRow<<<2, 2>>>(d_edge_up2, d_mat);

        matrixFill<<<1, 3>>>(d_seq1, d_seq2, d_mat, 4, 4);

        cudaMemcpy(mat, d_mat, 4*4*sizeof(int), cudaMemcpyDeviceToHost);
    }

    if (next == 1){
        initCol<<<2, 2>>>(d_edge_left2, d_mat);
        initRow<<<2, 2>>>(d_edge_up2+4, d_mat);

        matrixFill<<<1, 3>>>(d_seq1+4, d_seq2, d_mat, 4, 4);

        cudaMemcpy(mat, d_mat, 4*4*sizeof(int), cudaMemcpyDeviceToHost);
    }

    if (next == 2){
        initCol<<<2, 2>>>(d_edge_left2+4, d_mat);
        initRow<<<2, 2>>>(d_edge_up2, d_mat);

        matrixFill<<<1, 3>>>(d_seq1, d_seq2+4, d_mat, 4, 4);

        cudaMemcpy(mat, d_mat, 4*4*sizeof(int), cudaMemcpyDeviceToHost);
    }


   //initial cordinate
    x = 3;
    y = 3;

    index = 0;

    while (x!=0 and y!=0){
        if (mat[(y-1)*4+x] > mat[y*4+x-1]){
            temp_x = x-1;
            temp_y = y;
            temp_value = mat[y*4+x-1];
            temp_direction = 'U'; //up
        }
        else{
            temp_x = x;
            temp_y = y-1;
            temp_value = mat[(y-1)*4+x];
            temp_direction = 'L'; //left
        }
        if (mat[(y-1)*4+x-1] > temp_value){
            x = temp_x;
            y = temp_y;
            track[index] = temp_direction;
        }
        else{
            x = x-1;
            y = y-1;
            track[index] = 'D';
        }
        index = index+1;
    }
}


