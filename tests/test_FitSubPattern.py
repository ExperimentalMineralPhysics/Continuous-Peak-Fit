import unittest

from cpf import PeakFunctions as pf


class TestCreateNewparams(unittest.TestCase):
    def setUp(self):
        pass

    def test_CreateNewparams(self):
        num_peaks = 1
        dfour = [[1]]
        hfour = [[0.7]]
        wfour = [[0.6]]
        pfour = [[0.5]]
        bgfour = [[0.1, 0.1]]
        order_peak = []
        # order_peak = {}
        # order_peak['peak'] = {}
        # order_peak['peak'][0] = {}
        for i in range(num_peaks):
            order_peak.append({"symmetry": 1})
            # order_peak['peak'][i]['symmetry'] = 1.
        # print(order_peak[0]['symmetry'])
        testparam = pf.create_newparams(
            num_peaks, dfour, hfour, wfour, pfour, bgfour, order_peak
        )
        # print(testparam)
        assert testparam == {
            "background": [[0.1, 0.1]],
            "peak": [
                {
                    "d-space": 1,
                    "height": 0.7,
                    "width": 0.6,
                    "profile": 0.5,
                    "symmetry": 1,
                }
            ],
        }
