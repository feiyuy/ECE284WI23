#include <iostream>
#include <vector>
#include <fstream>
#include <string>
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

    std::fstream new_file;
    std::string sa;
    new_file.open(refFilename, std::ios::in); 
    getline(new_file, sa);
    char* s1 = sa.c_str();

    int i,l1;
    l1 = strlen(s1);
    int seq1[l1];
    for (i=0;i<l1;i++){
        switch(s1[i]) {
        case 'A':
            seq1[i] = 0;
            break;
        case 'C':
            seq1[i] = 1;
            break;
        case 'G':
            seq1[i] = 2;
            break;
        case 'T':
            seq1[i] = 3;
            break;
        default:
            seq1[i] = 0;
            break;
        }
    }

    getline(new_file, sa);
    char* s2 = sa.c_str();
    int j,l2;
    l2 = strlen(s2);
    int seq1[l1];
    for (j=0;i<l2;j++){
        switch(s2[j]) {
        case 'A':
            seq2[j] = 0;
            break;
        case 'C':
            seq2[j] = 1;
            break;
        case 'G':
            seq2[j] = 2;
            break;
        case 'T':
            seq2[j] = 3;
            break;
        default:
            seq2[j] = 0;
            break;
        }
    }

    new_file.close(); 

    int t,track;
    int dist[51][51];

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
            if(seq1[i-1] == seq2[j-1]) {
                track= 0;
            } else {
                track = 1;
            }
            t = MIN((dist[i-1][j]+1),(dist[i][j-1]+1));
            dist[i][j] = MIN(t,(dist[i-1][j-1]+track));
        }
    }
    std::cout<<"The Levinstein distance is:"<<dist[l2][l1];
    fprintf(stdout, "\nCompleted in %ld msec \n\n", timer.Stop());

    int mat[l1*l2];
    timer.Start();
    fprintf(stdout, "\nCompute Levinstein distance in GPU.\n");
    compute(seq1, seq2, l1, l2, mat);
    fprintf(stdout, "\nCompleted in %ld msec \n\n", timer.Stop());
    
    return 0;
}

