#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <boost/program_options.hpp>
#include <tbb/task_scheduler_init.h>
#include "seedTable.cuh"
#include "zlib.h"
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) < (y) ? (y) : (x))

// For parsing the command line values
namespace po = boost::program_options;

int main(int argc, char** argv) {
    // Timer below helps with the performance profiling (see timer.hpp for more
    // details)
    Timer timer;

    std::string refFilename;

    // Parse the command line options
    po::options_description desc{"Options"};
    desc.add_options()
    ("reference,r", po::value<std::string>(&refFilename)->required(), "Input FASTA file name [REQUIRED].")
    ("help,h", "Print help messages");

    po::options_description allOptions;
    allOptions.add(desc);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(allOptions).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << desc << std::endl;
        exit(1);
    }

    int i,j,l1,l2,t,track;
    int dist[51][51];
    //take the strings as input
    char s1[] = "GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGC";
    char s2[] = "CTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCT";
    //stores the lenght of strings s1 and s2
    l1 = strlen(s1);
    l2= strlen(s2);

    timer.Start();
    fprintf(stdout, "\nCompute Levinstein distance in CPU.\n");
    for(i=0;i<=l1;i++) {
        dist[0][i] = i;
    }
    for(j=0;j<=l2;j++) {
        dist[j][0] = j;
    }
    for (j=1;j<=l1;j++) {
        for(i=1;i<=l2;i++) {
            if(s1[i-1] == s2[j-1]) {
                track= 0;
            } else {
                track = 1;
            }
            t = MIN((dist[i-1][j]+1),(dist[i][j-1]+1));
            dist[i][j] = MIN(t,(dist[i-1][j-1]+track));
        }
    }
    std::cout<<"The Levinstein distance is:"<<dist[l2][l1];
    fprintf(stdout, "Completed in %ld msec \n\n", timer.Stop());

    int L1 = MIN(l1, l2);
    int L2 = MAX(l1, l2);

    int *d_mat;
    cudaMalloc(&d_mat, L1*L2*sizeof(int));

    initRows <<5, L1/5>>>(d_mat, L1);
    initCols <<5, L2/5>>>(d_mat, L1, L2);




    return 0;
}

