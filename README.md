# TS-EMO
This repository contains the source code for the “Thompson sampling efficient multiobjective optimization” (TSEMO) algorithm. 

To use TSEMO download all files contained in the repository and run the algorithm on the required test-function as shown in "TSEMO_Example".

______________________________________________________________________________________________________________________________
To cite TSEMO use the following publication (Bradford2018):

E. Bradford, A. M. Schweidtmann, and A. Lapkin, “Efficient multiobjective
optimization employing Gaussian processes, spectral sampling
and a genetic algorithm”, Journal of Global Optimization, pp. 1–33,
2018.

@article{Bradford2018,
author="Bradford, Eric
and Schweidtmann, Artur M.
and Lapkin, Alexei",
title="Efficient multiobjective optimization employing Gaussian processes, spectral sampling and a genetic algorithm",
journal="Journal of Global Optimization",
year="2018",
volume="71",
number="2",
pages="407--438",
doi="10.1007/s10898-018-0609-2"}

____________________________________________________________________________________________________________________________
The algorithm has already been successfully applied to several expensive multiobjective optimization problems:

- Determination of optimal conditions of a fully-automated chemical reactor system trading-off yield and environmental factors (Schweidtmann2018) 

- Optimization of a chemical process using a life-cycle assessment and cost simulation (Helmdach2017) 

@article{Schweidtmann2018,
  title={Machine learning meets continuous flow chemistry: Automated optimization towards the Pareto front of multiple objectives},
  author={Schweidtmann, Artur M and Clayton, Adam D and Holmes, Nicholas and Bradford, Eric and Bourne, Richard A and Lapkin, Alexei A},
  journal={Chemical Engineering Journal},
  year={2018},
  publisher={Elsevier}
}

@article{Helmdach2017,
  title={A Multiobjective Optimization Including Results of Life Cycle Assessment in Developing Biorenewables-Based Processes},
  author={Helmdach, Daniel and Yaseneva, Polina and Heer, Parminder K and Schweidtmann, Artur M and Lapkin, Alexei A},
  journal={ChemSusChem},
  volume={10},
  number={18},
  pages={3632--3643},
  year={2017},
  publisher={Wiley Online Library}
}
