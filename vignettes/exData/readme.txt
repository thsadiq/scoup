% ><>< ================================================================ ><>< %
% ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< %
% ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< %
% ><>< ================================================================ ><>< %

FIG1:
This folder contains a R script used to execute the code presented in the
first image. The times it took to complete 50 implementations of the code
are saved in the "times.csv" file. Note that the times are in seconds.

FIG2:
This folder contains the codon sequence data that were analysed to obtain the
outputs summarised in the second image of the article. The R script used to
execute the scoup simulation process is included in the folder, together with
a sample control file used to obtain the omega estimates in PAML 4.9. Many of
the input parameters of the simulator were fixed for all the data sets. These
include,
  - OU asymptotic mean = 0,
  - Branch length = 0.10 (for each of the branches on the balanced tree),
  - Terminal nodes = 8,
  - Number of codon sites = 1000,
  - Population size = 1000,
  - Sampling distribution for selection coefficients = Gamma(),
Four parameters were varied and these are listed below according to the file
names that contain the respective data.
  - vN1vS*ouS*ouT* = Variance of the amino acid coefficients input as zero,
  - vN2vS*ouS*ouT* = Variance of the amino acid coefficients input as 0.10,
  - vN*vS1ouS*ouT* = Variance of synonymous coefficients input as zero (i.e. vN/vS=0),
  - vN*vS2ouS*ouT* = Variance of synonymous coefficients input as 0.02 (i.e. vN/vS=5),
    + Note: Variance of synonymous coefficients input as 1/(vN/vS) when vN=0 (i.e. .2),
  - vN*vS*ouS1ouT* = OU asymptotic variance set as zero
  - vN*vS*ouS2ouT* = OU asymptotic variance set as 0.34
  - vN*vS*ouS3ouT* = OU asymptotic variance set as 0.67
  - vN*vS*ouS4ouT* = OU asymptotic variance set as 1.00
  - vN*vS*ouS*ouT1 = OU reversion parameter value fixed as 0.01
  - vN*vS*ouS*ouT2 = OU reversion parameter value fixed as 0.34
  - vN*vS*ouS*ouT3 = OU reversion parameter value fixed as 0.67
  - vN*vS*ouS*ouT4 = OU reversion parameter value fixed as 1.00
In addition, the folder contains a text file (evolTree.txt) where the phylogeny
used for the simulations is saved, two Excel files "dnds.csv" and "omega.csv"
that contain the analytical (dN/dS) and the inferred (ω) estimates of the
selection effect with respect to the simulated sequences.

% ><>< ================================================================ ><>< %
% ><>< ================================================================ ><>< %