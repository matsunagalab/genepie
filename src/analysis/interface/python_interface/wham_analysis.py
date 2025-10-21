import os
from . import genesis_exe


def test_wham_analysis():
    pmf = genesis_exe.wham_analysis(
            cvfile = '../../../../tests/regression_test/test_analysis/trajectories/triala_cv/{}.dis',
            dimension     = 1,
            nblocks       = 1,
            temperature   = 300.0,
            tolerance     = 10E-08,
            rest_function = (1, ),
            grids         = ((0.0, 15.0, 301), ),
            constant      =  ((1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2,
                               1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2,), ),
            reference     = ((1.80, 2.72, 3.64, 4.56, 5.48, 6.40, 7.32,
                              8.24, 9.16, 10.08, 11.00, 11.92, 12.84, 13.76,), ),
            is_periodic   = (False, ),
            )
    print(pmf)


def main():
    test_wham_analysis()


if __name__ == "__main__":
    main()
