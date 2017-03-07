import numpy as np
from credo.modelresult import ModelResult
from credo.systest import BaseWithinTolTC, calc_dist_errors
from xml.etree import ElementTree as etree

class DigitisedRadialFieldResult(ModelResult):
    """Digitised radial results for a field at a given time."""
    def __init__(self, modelName, fileName, field, outputIndex, ordering_map=None,
                 fieldname_map=None):
        from os.path import dirname
        super(DigitisedRadialFieldResult, self).__init__(modelName, dirname(fileName),
                                            ordering_map=ordering_map,
                                            fieldname_map=fieldname_map)
        self.field = field
        self.outputIndex = outputIndex
        self.data = np.loadtxt(fileName)
    def getRadii(self):
        return self.data[:, 0]
    def _getFieldAtOutputIndex(self, field, outputIndex):
        if field == self.field and outputIndex == self.outputIndex:
            return self.data[:, 1]
        else: return None
        
class RadialSolutionWithinTolTC(BaseWithinTolTC):
    """
    For testing against radial solution results. 
    """

    def __init__(self,
                 fieldsToTest = None,
                 defFieldTol = 0.01,
                 fieldTols = None,
                 expected = None,
                 absoluteErrorTol = 1.0,
                 testOutputIndex = -1,
                 maxRadius = None):
        BaseWithinTolTC.__init__(self,
                                 fieldsToTest = fieldsToTest,
                                 defFieldTol = defFieldTol,
                                 fieldTols = fieldTols,
                                 expected = expected,
                                 absoluteErrorTol = absoluteErrorTol)
        self.testOutputIndex = testOutputIndex
        self.maxRadius = maxRadius

    def _checkFieldWithinTol(self, field, mResult):
        fieldTol = self._getTolForField(field)
        r_result = np.array([pos[0] for pos in mResult.getPositions()])
        result = mResult.getFieldAtOutputIndex(field, self.testOutputIndex)
        if self.maxRadius:
            ir = np.where(r_result <= self.maxRadius)
            r_result = r_result[ir]
            result = result[ir]
        r_expected = self.expected.getRadii()
        expected = self.expected.getFieldAtOutputIndex(field, self.testOutputIndex)
        errors = calc_dist_errors(r_result, result,
                                  r_expected, expected,
                                  self.absoluteErrorTol, logx = True)
        fieldResult = all(e <= fieldTol for e in errors)
        return fieldResult, errors

    def _writeXMLCustomSpec(self, specNode):
        etree.SubElement(specNode, 'testOutputIndex',
            value = str(self.testOutputIndex))
        BaseWithinTolTC._writeXMLCustomSpec(self, specNode)
    

