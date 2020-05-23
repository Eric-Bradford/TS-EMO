# Thompson sampling efficient multiobjective optimization
This repository contains the source code for the “Thompson sampling efficient multiobjective optimization” (TSEMO) algorithm outlined in [(Bradford et al., 2018)](#Bradford2018). The algorithm is written to optimize expensive, black-box functions involving multiple conflicting criteria by employing Gaussian process surrogates. It is often able to determine a good approximation of the true Pareto front in signficantly less iterations than genetic algorithms. To cite TSEMO use [(Bradford et al., 2018)](#Bradford2018).

<img src="/Old_versions/Images/GP_sample_graphs.jpg" width="400">

## Getting started
To use TSEMO download all files contained in the repository and run the algorithm on the required test-function as shown in the example matlab file [TSEMO_Example](TSEMO_Example.m). To use the algorithm on your own functions simply copy the same format as the functions shown in the [test-function folder](/Test_functions/). The algorithm can be applied to any number of inputs and objectives. 

## Example applications
The algorithm has been successfully applied to several expensive multiobjective optimization problems:

* Determination of optimal conditions of a fully-automated chemical reactor system trading-off yield and environmental factors [(Schweidtmann et al., 2018)](#Schweidtmann2018) including multi-step reactions and separation processes [(Clayton et al., 2020)](#Clayton2020)

![](https://ars.els-cdn.com/content/image/1-s2.0-S1385894718312634-gr2.jpg)

* Optimization of a chemical process using a life-cycle assessment and cost simulation [(Helmdach et al., 2018)](#Helmdach2017) 

* Solvent selection for asymmetric catalysis using molecular descriptors [(Amar et al., 2019)](#Amar2019)


## References
E. Bradford, A. M. Schweidtmann, and A. A. Lapkin, [Efficient multiobjective optimization employing Gaussian processes, spectral sampling and a genetic algorithm](https://link.springer.com/article/10.1007/s10898-018-0609-2/), Journal of Global Optimization, vol. 71, no. 2, pp. 407–438, 2018.

<a name="Bradford2018">
</a>

A. M. Schweidtmann, A. D. Clayton, N. Holmes, E. Bradford, R. A. Bourne, and A. A. Lapkin, [Machine learning meets continuous flow chemistry: Automated optimization towards the Pareto front of multiple objectives](https://www.sciencedirect.com/science/article/pii/S1385894718312634), Chemical Engineering Journal, vol. 352, pp. 277-282, 2018.    

<a name="Schweidtmann2018">
</a>

D. Helmdach, P. Yaseneva, K. P. Heer, A. M. Schweidtmann, and A. A. Lapkin, [A Multiobjective Optimization Including Results of Life Cycle Assessment in Developing Biorenewables-Based Processes](https://onlinelibrary.wiley.com/doi/abs/10.1002/cssc.201700927), ChemSusChem, vol. 10, no. 18, pp. 3632-3643, 2017.  

<a name="Helmdach2017">
</a>

Y. Amar, A. M. Schweidtmann, P. Deutsch, L. Cao, and A. A. Lapkin, [Machine learning and molecular descriptors enable rational solvent selection in asymmetric catalysis](https://pubs.rsc.org/en/content/articlelanding/2019/sc/c9sc01844a#!divAbstract), Chemical Science, vol. 10, no. 27, pp. 6697-6706, 2019. 
<a name="Amar2019">
</a>

A. Clayton, A. M. Schweidtmann, G. Clemens, J. Manson, C. Taylor, C. Nino, T. Chamberlain, N. Kapur, A. Blacker, A. A. Lapkin, R. Bourne [Automated self-optimisation of multi-step reaction and separation processes using machine learning](https://doi.org/10.1016/j.cej.2019.123340), Chemical Engineering Journal, vol. 384, 123340, 2020. 
<a name="Clayton2020">
</a>
