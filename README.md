FiniteInflation
===============

Code for Cosmoslik with Class for MCMC

This code runs a modified power spectrum using a kinetic dominated start for inflation. See http://arxiv.org/abs/1503.04872 for more details.

Class was downloaded from: https://github.com/lesgourg/class_public/releases/tag/v2.4.1
Cosmoslik is from: https://github.com/marius311/cosmoslik


Setup instructions
==================

Git clone this repo, then run from the CLASS directory:

make clean
make all 

add to .bashrc the locations of likelihood files and also the main cosmoslik directory so python knows where to look:

export PYTHONPATH=/.../cosmomc/likelihoods/clik_0313/plc-1.0/lib/python2.7/site-packages:$PYTHONPATH
export PYTHONPATH=~/.../cosmoslik:$PYTHONPATH

add to .bash_profile:

if [ -f ~/.bashrc ]; then
   source ~/.bashrc
fi

Configure and build python code (run from cosmoslik directory):

./waf configure --install_all_deps

./waf install

Add to .bash_profile:

if [ -f ~/planck_data/plc-1.0/bin/clik_profile.sh ]; then
source ~/planck_data/plc-1.0/bin/clik_profile.sh
fi



To run the chain:

nohup mpiexec -n 5 python -m cosmoslik myprogram.py --traceback &

where 5 = 4 + 1 corresponds to 4 chains and one master process, and myprogram.py is replaced with the filename of the parameter settings code.
