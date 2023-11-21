## Code for the paper 'Approximate Bayesian inference from noisy likelihoods with Gaussian process emulated MCMC'

This repository contains a MATLAB implementation of the GP-MH (and MH-BLFI) algorithm used in the paper <https://arxiv.org/abs/2104.03942>. Also scripts for drawing all the illustrative figures in the paper are provided. The implementation contains also some extra features not used in the paper (e.g. additional test problems and alternative acquisition functions) which are, however, not carefully tested.

The code is partially based on the code of our earlier paper 'Parallel Gaussian process surrogate Bayesian inference with noisy likelihood evaluations' in <https://github.com/mjarvenpaa/parallel-GP-SL>.

Note that the code implementation is provided mainly to demonstrate the promise of the algorithm and for reproducibility and should not be considered as mature, easy-to-use inference software in its current version.

## Getting started

Check out and run the file `GPMH_test_run`. The files in the `extra` folder starting with `demo_` can be used to draw the figures of the paper and additional illustrations.

## Installation and external code

Place the code files to your working directory and make sure that all the folders are contained in the MATLAB search path. The list below contains links to some external software/code needed. Just obtain these code packages and place them to the working directory.

* GPstuff (used e.g. for MAP estimation of the GP hyperparameters): <https://github.com/gpstuff-dev/gpstuff>
* DRAM MCMC (used e.g. for sampling in MH-BLFI): <http://helios.fmi.fi/~lainema/dram/>
* Cell biology and g-and-k models: <https://github.com/cdrovandi/Bayesian-Synthetic-Likelihood>
* `export_fig` (needed for saving some figures nicely): <https://github.com/altmany/export_fig>
* `subtightplot` (needed for plotting some figures nicely): <https://se.mathworks.com/matlabcentral/fileexchange/39664-subtightplot>
* `shadedplot` (needed for plotting some figures nicely): <https://se.mathworks.com/matlabcentral/fileexchange/18738-shaded-area-plot>
* `shadederrorbar` (needed for plotting some figures nicely): <https://se.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar>

The code for bacterial infections and Lorenz models is also not included here but can be obtained from the author on request.

Important: DRAM package contains a file `mad.m` which should not shadow the MATLAB file with the same name, otherwise an error occurs. One thus needs to rename (or remove) `mad.m` contained in the DRAM or ensure otherwise that MATLABs `mad.m` is always called.

If you want to try DIRECT or CMAES algorithms for optimising the acquisition functions you can obtain implementations from <https://github.com/npinto/direct> and <http://cma.gforge.inria.fr/> or from the author on request.

## Ideas for possible extensions/improvements

* Better modelling of the variance of the noisy log-likelihood evaluations.
* The accuracy of the MH accept/reject decision is currently controlled via an upper bound $\varepsilon$. Alternative strategies that e.g. take into account the distance between the proposal and the current point of the MH chain or are based on some notion of average accuracy could be investigated.
* The current GP-MH implementation might still sometimes get stuck or fail otherwise due to poor GP fits. Robustness could probably be improved by better engineering or theoretical considerations.
* Parallel acquisitions.

## Support

If you have questions or would like to use the provided implementation or methodology for your own inference problem, please contact <m.j.jarvenpaa@medisin.uio.no> or <jarvenpaamj@gmail.com>.

## License

This software is distributed under the GNU General Public Licence (version 3 or later); please refer to the file LICENSE.txt for details.
