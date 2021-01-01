import numpy as np
from AbundanceMatching import AbundanceFunction, LF_SCATTER_MULT

_lf_table = '''
-24.7 -6.285
-24.5 -5.861
-24.3 -5.518
-24.1 -5.161
-23.9 -4.903
-23.7 -4.651
-23.5 -4.411
-23.3 -4.170
-23.1 -3.953
-22.9 -3.751
-22.7 -3.555
-22.5 -3.383
-22.3 -3.229
-22.1 -3.095
-21.9 -2.971
-21.7 -2.876
-21.5 -2.801
-21.3 -2.722
-21.1 -2.668
-20.9 -2.604
-20.7 -2.559
-20.5 -2.509
-20.3 -2.494
-20.1 -2.487
-19.9 -2.477
-19.7 -2.451
-19.5 -2.447
-19.3 -2.409
-19.1 -2.399
-18.9 -2.377
-18.7 -2.371
-18.5 -2.348
-18.3 -2.305
-18.1 -2.309
-17.9 -2.335
-17.7 -2.350'''.strip().replace('\n', ' ')

_lf = np.fromstring(_lf_table, sep=' ').reshape(-1, 2)

_lf[:, 1] = 10.0**_lf[:, 1]


def test_AbundanceFunction():
    AbundanceFunction(*_lf.T, ext_range=(-25, -16))


def test_deconv():
    af = AbundanceFunction(*_lf.T, ext_range=(-25, -16))
    af.deconvolute(0.2 * LF_SCATTER_MULT, 20)


if __name__ == "__main__":
    test_AbundanceFunction()
    test_deconv()
