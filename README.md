# Sarcomeric_pattern_formation
This GitHub repository includes codes and examples for the image-analysis algorithm and simulation code used in the PhD thesis "Physical principles of pattern formation during myofibrillogenesis" by Francine Kolley, TU Dresden, and the corresponding manuscript (https://journals.aps.org/prxlife/abstract/10.1103/PRXLife.2.013002). If you use the code, please cite the manuscript.

This repository contains the following supporting code:

(1) A tracking-free algorithm to compute correlation functions from multi-channel images, using a steerable filter to compute the local nematic director.
    The algorithm was modified in collaboration with Benoit Dehapiot based on an older version that can be found at https://github.com/BDehapiot/BDProject_FDomMuscle_2D .
    The code uses a third-party package for the computatio of the steerable filter for which a separate licence applies. 
    This repository includes a verbatim copy of this package, which was downloaded from the following URL http://www.francoisaguet.net/software.html .
    
(2) Example code for our minimal model I describing sarcomeric pattern formation driven by non-local interactions. 
    Code examples are given for both mean-field simulations and agent-based simulations. 
