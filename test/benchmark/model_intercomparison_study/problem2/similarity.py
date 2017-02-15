import numpy as np
from credo.modelresult import ModelResult
from credo.systest import BaseWithinTolTC, calc_dist_errors
from xml.etree import ElementTree as etree

class DigitisedSimilarityResult(ModelResult):
    """Digitised similarity solution results for a set of fields."""
    def __init__(self, modelName, fileNames, ordering_map = None,
                 fieldname_map = None):
        super(DigitisedSimilarityResult, self).__init__(modelName, '',
                                            ordering_map = ordering_map,
                                            fieldname_map = fieldname_map)
        self.data = {}
        for field, fileName in fileNames.items():
            self.data[field] = np.loadtxt(fileName)
    def getSimilarityVariables(self, field): return self.data[field][:,0]
    def getSimilarityValues(self, field): return self.data[field][:,1]

class SimilaritySolutionWithinTolTC(BaseWithinTolTC):
    """
    For testing against semi-analytical similarity solution results. 
    """

    def __init__(self,
                 fieldsToTest = None,
                 defFieldTol = 0.01,
                 fieldTols = None,
                 expected = None,
                 absoluteErrorTol = 1.0,
                 testCellIndices = [0]):
        BaseWithinTolTC.__init__(self,
                                 fieldsToTest = fieldsToTest,
                                 defFieldTol = defFieldTol,
                                 fieldTols = fieldTols,
                                 expected = expected,
                                 absoluteErrorTol = absoluteErrorTol)
        self.testCellIndices = testCellIndices

    def _checkFieldWithinTol(self, field, mResult):

        fieldTol = self._getTolForField(field)
        def sort_arrays(sims, results):
            results = np.array(results)
            sims = np.array(sims)
            isort = np.argsort(sims)
            return sims[isort], results[isort]

        results, result_sims = [], []
        for testCellIndex in self.testCellIndices:
            result = mResult.getFieldHistoryAtCell(field, testCellIndex)
            results += list(result)
            result_times = mResult.getTimes()
            pos = mResult.getPositions()[testCellIndex]
            r2 = pos[0] * pos[0]
            result_sim = result_times / r2
            result_sims += list(result_sim)
        result_sims, results = sort_arrays(result_sims, results)

        if hasattr(self.expected, 'getSimilarityVariables'):
            expected_sims = self.expected.getSimilarityVariables(field)
            expecteds = self.expected.getSimilarityValues(field)
        else:
            expecteds, expected_sims = [], []
            for testCellIndex in self.testCellIndices:
                pos = mResult.getPositions()[testCellIndex]
                if callable(self.expected):
                    expected = np.array([self.expected(pos, t)
                                         for t in result_times])
                    expected_times = result_times
                else:
                    expected = self.expected.getFieldHistoryAtCell(field,
                                                                   testCellIndex)
                    expected_times = self.expected.getTimes()
                expecteds += list(expected)
                r2 = pos[0] * pos[0]
                expected_sim = expected_times / r2
                expected_sims += list(expected_sim)
            expected_sims, expecteds = sort_arrays(expected_sims,
                                                       expecteds)

        errors = calc_dist_errors(result_sims, results,
                                  expected_sims, expecteds,
                                  self.absoluteErrorTol, logx = True)

        fieldResult = all(e <= fieldTol for e in errors)
        return fieldResult, errors

    def _writeXMLCustomSpec(self, specNode):
        etree.SubElement(specNode, 'testCellIndices',
            value = str(self.testCellIndices))
        BaseWithinTolTC._writeXMLCustomSpec(self, specNode)
    
