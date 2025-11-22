from paraview.util.vtkAlgorithm import *

@smproxy.filter()
@smproperty.input(name="Input")
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=True)
class SimplePassThrough(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
            nInputPorts=1, nOutputPorts=1, inputType='vtkDataSet')

    def RequestData(self, request, inInfo, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkDataSet
        input = vtkDataSet.GetData(inInfo[0], 0)
        output = vtkDataSet.GetData(outInfo, 0)
        output.ShallowCopy(input)
        return 1
