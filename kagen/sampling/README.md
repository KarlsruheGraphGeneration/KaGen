## Efficient Random Sampling  [![Travis-CI Status](https://travis-ci.org/lorenzhs/sampling.svg?branch=master)](https://travis-ci.org/lorenzhs/sampling)
#### Parallel, Vectorized, Cache-Efficient, and Online

This is the code to accompany our eponymous paper: *Sanders, P., Lamm, S., Hübschle-Schneider, L., Schrade, E. and Dachsbacher, C., 2018. Efficient Parallel Random Sampling—Vectorized, Cache-Efficient, and Online. ACM Transactions on Mathematical Software (TOMS), 44.3 (2018): 29.*

You can find the freely accessible author's version [in the arXiv](https://arxiv.org/abs/1610.05141).

If you use this library in the context of an academic publication, we ask that you cite our paper:
```bibtex
@article{sanders2018sampling,
  title = {Efficient Parallel Random Sampling---Vectorized, Cache-Efficient, and Online},
  author = {Sanders, Peter and Lamm, Sebastian and H{\"u}bschle-Schneider, Lorenz
            and Schrade, Emanuel and Dachsbacher, Carsten},
  journal = {ACM Trans. Math. Softw.},
  year = {2018},
  issn = {0098-3500},
  volume = {44},
  number = {3},
  articleno = {29},
  pages = {29:1--29:14},
  issue_date = {April 2018},
  publisher = {ACM},
  doi = {10.1145/3157734}
}
```
### Building

Build with cmake (version 2.8.12 or later is required). Remember to fetch the submodules before compiling: `git submodule update --init`. A compiler compatible with C++11 is required.

**Optional Dependencies:** An MPI library for Algorithm P's test runner (Algorithm P itself is implemented independently of MPI), Intel Math Kernel Library (MKL) for faster random variate generation

### Tests

To run the tests, set the cmake variable `SAMPLING_BUILD_TESTS` to `ON`.  Make sure that the `extlib/googletest` submodule is present.  Build, then execute the tests with `ctest`.

We plan to make the test suite more comprehensive in the future.

**[License](/LICENSE):** 2-clause BSD
