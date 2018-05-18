# Probabilistic Approaches to the AXB = YCZ Calibration Problem in Multi-Robot Systems
This article introduces two probabilistic solvers for the AXB=YCZ calibration problem in multi-robot systems, where A, B, C are time-varying rigid body transformations measured from sensors and X, Y, Z are unknown static transformations to be calibrated. Comparisons with other solvers have been made and the proposed probabilistic solvers are outperformed, especially in the case with the lack of correspondance of the measured data. The robustness of the proposed solvers is further improved by hybrid method and iterative refinement. This repository contains the Matlab implementations for probabilistic approaches to the AXB=YCZ calibrations problem. The work has been published in the journal of Autonomous Robots [[Link to the article]](https://rpk.lcsr.jhu.edu/publications/#Sensor_Calibration).

Authors: Qianli Ma (<mqianli1@jhu.edu>), Zachariah Goh (<zach_goh@yahoo.com>), Sipu Ruan (<ruansp@jhu.edu>), Gregory S. Chirikjian (<gchirik1@jhu.edu>).

## Repository Structure
The codes in this repository are organized to be self-contained, with descriptions as follows:

1. Test files are located in the root directory, with prefix "main_": \
  (1) "main_comparison_cov.m": solvers comparisons for synthetic data, with the change of covariance of noise; \
  (2) "main_comparison_disorder.m": comparisons for synthetic data that is scrambled; \
  (3) "main_NAO_data_analysis.m": comparions for real data collected from two NAO humanoid robot.

2. Functions that implements different solvers are under the folder "solvers/": \
  (1) "axbyczProb1.m": proposed probabilistic solver 1; \
  (2) "axbyczProb2.m": proposed probabilistic solver 2; \
  (3) "axbyczProb3.m": Robust iterative refinement for the proposed probabilistic solvers; \
  (4) Other functions: implementations of other solvers for comparisons.

3. Real experimental data collected from two humanoid robots are gathered in the folder "real_data/", with one sample trial: \
  (1) "transform_ABC_unified_fixA.mat": measured data for rigid body transforms B, C when A is fixed; \
  (2) "transform_ABC_unified_fixC.mat": measured data for rigid body transforms A, B when C is fixed; \
  (3) "transform_ABC_unified.mat": combined data of the previous two sets.

4. Supporting functions for running the solvers and test files are under "util/" folder.

## Citations
To reference our work in publications, please cite:

<cite>Ma, Q., Goh, Z., Ruan, S. and Chirikjian, G.S., 2018. Probabilistic approaches to the AXB= YCZ AXB= YCZ calibration problem in multi-robot systems. Autonomous Robots, pp.1-24.</cite>

```
@article{ma2018probabilistic,
  title={Probabilistic approaches to the AXB=YCZ calibration problem in multi-robot systems},
  author={Ma, Qianli and Goh, Zachariah and Ruan, Sipu and Chirikjian, Gregory S},
  journal={Autonomous Robots},
  pages={1--24},
  year={2018},
  publisher={Springer}
}
```
