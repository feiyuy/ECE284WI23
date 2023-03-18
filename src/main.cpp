#include <iostream>
#include <vector>
#include <boost/program_options.hpp>
#include <tbb/task_scheduler_init.h>
#include "seedTable.cuh"
#include "kseq.h"
#include "zlib.h"

// For parsing the command line values
namespace po = boost::program_options;

// For reading in the FASTA file
KSEQ_INIT2(, gzFile, gzread)

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

    // Read input sequence as kseq_t object
    timer.Start();
    fprintf(stdout, "Reading input sequence and compressing to two-bit encoding.\n");
    gzFile fp = gzopen(refFilename.c_str(), "r");
    if (!fp) {
        fprintf(stdout, "ERROR: Cannot open file: %s\n", refFilename.c_str());
        exit(1);
    }
    kseq_t *record = kseq_init(fp);
    int n;
    if ((n = kseq_read(record)) < 0) {
        fprintf(stdout, "ERROR: No records found!\n");
        exit(1);
    }
    printf("Sequence name: %s\n", record->name.s);
    printf("Sequence size: %zu\n", record->seq.l);
    fprintf(stdout, "Completed in %ld msec \n\n", timer.Stop());


    return 0;
}

