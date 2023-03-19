# ECE 284 Final Project

## Contents
- [ECE 284 Final Project](#ece-284-final-project)
  - [Contents](#contents)
  - [Project Topic](#project-topic)
  - [Setting Up](#setting-up)
  - [Running the Code](#running-the-code)

## Project Topic

Using GPUs to Speed-Up Levenshtein Edit Distance Computation


## Setting Up

1. SSH into the DSMLP server (dsmlp-login.ucsd.edu) using the AD account. MacOS and Linux users can SSH into the server using the following command:

```
ssh $username@dsmlp-login.ucsd.edu
```

2. Clone the assignment repository into $HOME directory using the following command:
```
cd ~
git clone https://github.com/feiyuy/ECE284WI23.git
```

## Running the Code

1. In this project we use file under data folder as input. By default the program will take the first 2 lines of the file as reference and target sequence to compute the Levinstein distance. Changes could be made to the input file to change the sequence length. Timer function is called to measure time comsuption. There is also a demo function that will take two sequences with length of 8 to do the backtracking using the segemented way to save memory.
 
2. We will be using a Docker container, namely `yatisht/ece284-wi23:latest`, for submitting a job on the cluster containing the right virtual environment to build and test the code. This container already contains the correct Cmake version, CUDA and Boost libraries preinstalled within Ubuntu-20.04 OS. Note that these Docker containers use the same filesystem as the DSMLP login server, and hence the files written to or modified by the conainer is visiable to the login server and vice versa. To submit a job that executes `run-commands.sh` script located inside the `ECE284WI23` direcotry on a VM instance with 8 CPU cores, 16 GB RAM and 1 GPU device (this is the maxmimum allowed request on the DSMLP platform), the following command can be executed from the VS Code or DSMLP Shell Terminal (replace the username and directory names below appropriately):

```
ssh $username@dsmlp-login.ucsd.edu /opt/launch-sh/bin/launch.sh -c 8 -g 1 -m 16 -i yatisht/ece284-wi23:latest -f ${HOME}/ECE284WI23/run-commands.sh
```

