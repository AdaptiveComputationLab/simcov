# state file generated using paraview version 5.7.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

##### INSERT BEGIN
import glob
import os
import sys
import argparse
import functools

argparser = argparse.ArgumentParser()
argparser.add_argument("--data", required=True, help='Specify the output directory containing the time series data')
options = argparser.parse_args()

def get_sample_fnames(prefix):
    global options
    fnames = sorted(glob.glob(options.data + '/' + prefix + '*'), key=lambda s: int(s.split('_')[2].split('.')[0]))
    return fnames


xdim = 0
ydim = 0
with open(options.data + '/sample_epicell_0.vtk') as f:
    for line in f.readlines():
        if line.startswith('DIMENSIONS'):
            xdim, ydim, _ = line.split()[1:]
            xdim = int(xdim) - 1
            ydim = int(ydim) - 1
            print('Dimensions of simulation are:', xdim, ydim)
            break
    else:
        print('Cannot find dimensions in sample files', file=sys.stderr)
        os.abort()

##### INSERT END

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Quartile Chart View'
cytokinechart = CreateView('QuartileChartView')
cytokinechart.ViewSize = [542, 278]
cytokinechart.LegendPosition = [260, 167]
cytokinechart.LeftAxisRangeMaximum = 0.05
cytokinechart.BottomAxisRangeMaximum = 200.0
cytokinechart.RightAxisRangeMaximum = 6.66
cytokinechart.TopAxisRangeMaximum = 6.66

# Create a new 'Quartile Chart View'
epicellchart = CreateView('QuartileChartView')
epicellchart.ViewSize = [541, 281]
epicellchart.ChartTitle = 'epicell averages'
epicellchart.ShowLegend = 0
epicellchart.LegendPosition = [214, 149]
epicellchart.LeftAxisRangeMaximum = 0.025
epicellchart.BottomAxisRangeMaximum = 200.0
epicellchart.RightAxisRangeMaximum = 6.66
epicellchart.TopAxisRangeMaximum = 6.66

# Create a new 'Render View'
epicellview = CreateView('RenderView', 'epicell view')
epicellview.ViewSize = [1093, 662]
epicellview.InteractionMode = 'Selection'
epicellview.AxesGrid = 'GridAxes3DActor'
epicellview.OrientationAxesVisibility = 0
##### REPLACE BEGIN
epicellview.CenterOfRotation = [2.0 * xdim, 2.0 * ydim, 0.5]
epicellview.StereoType = 'Crystal Eyes'
epicellview.CameraPosition = [2.0 * xdim, 2.0 * ydim, 14 * max(xdim, ydim)]
epicellview.CameraFocalPoint = [2.0 * xdim, 2.0 * ydim, 0.5]
##### REPLACE END
epicellview.CameraFocalDisk = 1.0
epicellview.CameraParallelScale = 42.875511494162964
epicellview.Background = [0.32, 0.34, 0.43]

# init the 'GridAxes3DActor' selected for 'AxesGrid'
epicellview.AxesGrid.Visibility = 1

# Create a new 'Quartile Chart View'
tcellchart = CreateView('QuartileChartView')
tcellchart.ViewSize = [541, 281]
tcellchart.ChartTitle = 't-cell average'
tcellchart.ShowLegend = 0
tcellchart.LegendPosition = [261, 181]
tcellchart.LeftAxisRangeMaximum = 0.007
tcellchart.BottomAxisRangeMaximum = 200.0
tcellchart.RightAxisRangeMaximum = 6.66
tcellchart.TopAxisRangeMaximum = 6.66

# Create a new 'Render View'
tcellview = CreateView('RenderView')
tcellview.ViewSize = [1092, 659]
tcellview.InteractionMode = 'Selection'
tcellview.AxesGrid = 'GridAxes3DActor'
tcellview.OrientationAxesVisibility = 0
##### REPLACE BEGIN
tcellview.CenterOfRotation = [2.0 * xdim, 2.0 * ydim, 0.5]
tcellview.StereoType = 'Crystal Eyes'
tcellview.CameraPosition = [2.0 * xdim, 2.0 * ydim, 14 * max(xdim, ydim)]
tcellview.CameraFocalPoint = [2.0 * xdim, 2.0 * ydim, 0.5]
##### REPLACE END
tcellview.CameraFocalDisk = 1.0
tcellview.CameraParallelScale = 42.875511494162964
tcellview.Background = [0.32, 0.34, 0.43]

# init the 'GridAxes3DActor' selected for 'AxesGrid'
tcellview.AxesGrid.Visibility = 1

# Create a new 'Quartile Chart View'
viruschart = CreateView('QuartileChartView')
viruschart.ViewSize = [541, 278]
viruschart.ChartTitle = 'virus average'
viruschart.ShowLegend = 0
viruschart.LegendPosition = [299, 184]
viruschart.LeftAxisRangeMaximum = 0.02
viruschart.BottomAxisRangeMaximum = 200.0
viruschart.RightAxisRangeMaximum = 6.66
viruschart.TopAxisRangeMaximum = 6.66

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.SplitVertical(1, 0.692385)
layout1.AssignView(3, epicellview)
layout1.SplitHorizontal(4, 0.500000)
layout1.AssignView(9, cytokinechart)
layout1.AssignView(10, viruschart)
layout1.SplitVertical(2, 0.689379)
layout1.AssignView(5, tcellview)
layout1.SplitHorizontal(6, 0.500000)
layout1.AssignView(13, tcellchart)
layout1.AssignView(14, epicellchart)

# ----------------------------------------------------------------
# restore active view
SetActiveView(epicellchart)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
legacyVTKReader2 = LegacyVTKReader(FileNames=get_sample_fnames('sample_chemokine_'))

# create a new 'Calculator'
calculator1 = Calculator(Input=legacyVTKReader2)
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'chemokine'
calculator1.Function = 'chemokine/255'

# create a new 'Plot Data Over Time'
plotDataOverTime1 = PlotDataOverTime(Input=calculator1)
plotDataOverTime1.FieldAssociation = 'Cells'

# create a new 'Legacy VTK Reader'
legacyVTKReader1 = LegacyVTKReader(FileNames=get_sample_fnames('sample_epicell_'))

# create a new 'Calculator'
calculator3 = Calculator(Input=legacyVTKReader1)
calculator3.AttributeType = 'Cell Data'
calculator3.ResultArrayName = 'dead epicells'
calculator3.Function = 'floor(epicell/4)'

# create a new 'Legacy VTK Reader'
legacyVTKReader3 = LegacyVTKReader(FileNames=get_sample_fnames('sample_virus_'))

# create a new 'Calculator'
calculator2 = Calculator(Input=legacyVTKReader3)
calculator2.AttributeType = 'Cell Data'
calculator2.ResultArrayName = 'virus'
calculator2.Function = 'virus/255'

# create a new 'Plot Data Over Time'
plotDataOverTime2 = PlotDataOverTime(Input=calculator2)
plotDataOverTime2.FieldAssociation = 'Cells'

# create a new 'Legacy VTK Reader'
legacyVTKReader4 = LegacyVTKReader(FileNames=get_sample_fnames('sample_tcelltissue_'))

# create a new 'Legacy VTK Reader'
legacyVTKReader5 = LegacyVTKReader(FileNames=get_sample_fnames('sample_icytokine_'))

# create a new 'Calculator'
calculator7 = Calculator(Input=legacyVTKReader5)
calculator7.AttributeType = 'Cell Data'
calculator7.ResultArrayName = 'icytokine'
calculator7.Function = 'icytokine/255'

# create a new 'Plot Data Over Time'
plotDataOverTime5 = PlotDataOverTime(Input=calculator7)
plotDataOverTime5.FieldAssociation = 'Cells'

# create a new 'Calculator'
calculator4 = Calculator(Input=calculator3)
calculator4.AttributeType = 'Cell Data'
calculator4.ResultArrayName = 'incubating epicells'
calculator4.Function = 'floor(2/(epicell+1))*epicell'

# create a new 'Calculator'
calculator5 = Calculator(Input=calculator4)
calculator5.AttributeType = 'Cell Data'
calculator5.ResultArrayName = 'expressing epicells'
calculator5.Function = 'floor(epicell/2*floor(3/(epicell+1)))'

# create a new 'Calculator'
calculator6 = Calculator(Input=calculator5)
calculator6.AttributeType = 'Cell Data'
calculator6.ResultArrayName = 'apoptotic epicells'
calculator6.Function = 'floor(epicell/3)*floor(4/(epicell+1))'

# create a new 'Plot Data Over Time'
plotDataOverTime4 = PlotDataOverTime(Input=legacyVTKReader4)
plotDataOverTime4.FieldAssociation = 'Cells'

# create a new 'Plot Data Over Time'
plotDataOverTime3 = PlotDataOverTime(Input=calculator6)
plotDataOverTime3.FieldAssociation = 'Cells'

# ----------------------------------------------------------------
# setup the visualization in view 'cytokinechart'
# ----------------------------------------------------------------

# show data from plotDataOverTime5
plotDataOverTime5Display = Show(plotDataOverTime5, cytokinechart)

# trace defaults for the display properties.
plotDataOverTime5Display.AttributeType = 'Row Data'
plotDataOverTime5Display.UseIndexForXAxis = 0
plotDataOverTime5Display.XArrayName = 'Time'
plotDataOverTime5Display.SeriesVisibility = ['icytokine (stats)', 'icytokine norm (stats)']
plotDataOverTime5Display.SeriesLabel = ['icytokine (stats)', 'icytokine (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)', 'icytokine norm (stats)', 'icytokine norm (stats)']
plotDataOverTime5Display.SeriesColor = ['icytokine (stats)', '0.666667', '0.333333', '1', 'N (stats)', '0.220005', '0.489998', '0.719997', 'Time (stats)', '0.300008', '0.689998', '0.289998', 'vtkValidPointMask (stats)', '0.6', '0.310002', '0.639994', 'icytokine norm (stats)', '0.333333', '0', '1']
plotDataOverTime5Display.SeriesPlotCorner = ['N (stats)', '0', 'Time (stats)', '0', 'icytokine (stats)', '0', 'icytokine norm (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime5Display.SeriesLabelPrefix = ''
plotDataOverTime5Display.SeriesLineStyle = ['N (stats)', '1', 'Time (stats)', '1', 'icytokine (stats)', '1', 'icytokine norm (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotDataOverTime5Display.SeriesLineThickness = ['N (stats)', '2', 'Time (stats)', '2', 'icytokine (stats)', '2', 'icytokine norm (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotDataOverTime5Display.SeriesMarkerStyle = ['N (stats)', '0', 'Time (stats)', '0', 'icytokine (stats)', '0', 'icytokine norm (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime5Display.ShowQuartiles = 0
plotDataOverTime5Display.ShowRanges = 0

# show data from plotDataOverTime1
plotDataOverTime1Display = Show(plotDataOverTime1, cytokinechart)

# trace defaults for the display properties.
plotDataOverTime1Display.AttributeType = 'Row Data'
plotDataOverTime1Display.UseIndexForXAxis = 0
plotDataOverTime1Display.XArrayName = 'Time'
plotDataOverTime1Display.SeriesVisibility = ['chemokine (stats)', 'chemokine norm (stats)']
plotDataOverTime1Display.SeriesLabel = ['chemokine (stats)', 'chemokine (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)', 'chemokine norm (stats)', 'chemokine norm (stats)']
plotDataOverTime1Display.SeriesColor = ['chemokine (stats)', '0', '0.666667', '1', 'N (stats)', '0.220005', '0.489998', '0.719997', 'Time (stats)', '0.300008', '0.689998', '0.289998', 'vtkValidPointMask (stats)', '0.6', '0.310002', '0.639994', 'chemokine norm (stats)', '1', '0.666667', '0']
plotDataOverTime1Display.SeriesPlotCorner = ['N (stats)', '0', 'Time (stats)', '0', 'chemokine (stats)', '0', 'chemokine norm (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime1Display.SeriesLabelPrefix = ''
plotDataOverTime1Display.SeriesLineStyle = ['N (stats)', '1', 'Time (stats)', '1', 'chemokine (stats)', '1', 'chemokine norm (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotDataOverTime1Display.SeriesLineThickness = ['N (stats)', '2', 'Time (stats)', '2', 'chemokine (stats)', '2', 'chemokine norm (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotDataOverTime1Display.SeriesMarkerStyle = ['N (stats)', '0', 'Time (stats)', '0', 'chemokine (stats)', '0', 'chemokine norm (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime1Display.ShowQuartiles = 0
plotDataOverTime1Display.ShowRanges = 0

# ----------------------------------------------------------------
# setup the visualization in view 'epicellchart'
# ----------------------------------------------------------------

# show data from plotDataOverTime3
plotDataOverTime3Display = Show(plotDataOverTime3, epicellchart)

# trace defaults for the display properties.
plotDataOverTime3Display.AttributeType = 'Row Data'
plotDataOverTime3Display.UseIndexForXAxis = 0
plotDataOverTime3Display.XArrayName = 'Time'
plotDataOverTime3Display.SeriesVisibility = ['apoptotic epicells (stats)', 'dead epicells (stats)', 'expressing epicells (stats)', 'incubating epicells (stats)']
plotDataOverTime3Display.SeriesLabel = ['apoptotic epicells (stats)', 'apoptotic epicells (stats)', 'dead epicells (stats)', 'dead epicells (stats)', 'epicell (stats)', 'epicell (stats)', 'expressing epicells (stats)', 'expressing epicells (stats)', 'incubating epicells (stats)', 'incubating epicells (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)']
plotDataOverTime3Display.SeriesColor = ['apoptotic epicells (stats)', '0.533333', '0', '0', 'dead epicells (stats)', '0', '0', '0', 'epicell (stats)', '0.220005', '0.489998', '0.719997', 'expressing epicells (stats)', '1', '0.333333', '0', 'incubating epicells (stats)', '0.666667', '0.666667', '0.498039', 'N (stats)', '1', '0.500008', '0', 'Time (stats)', '0.650004', '0.340002', '0.160006', 'vtkValidPointMask (stats)', '0', '0', '0']
plotDataOverTime3Display.SeriesPlotCorner = ['N (stats)', '0', 'Time (stats)', '0', 'apoptotic epicells (stats)', '0', 'dead epicells (stats)', '0', 'epicell (stats)', '0', 'expressing epicells (stats)', '0', 'incubating epicells (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime3Display.SeriesLabelPrefix = ''
plotDataOverTime3Display.SeriesLineStyle = ['N (stats)', '1', 'Time (stats)', '1', 'apoptotic epicells (stats)', '1', 'dead epicells (stats)', '1', 'epicell (stats)', '1', 'expressing epicells (stats)', '1', 'incubating epicells (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotDataOverTime3Display.SeriesLineThickness = ['N (stats)', '2', 'Time (stats)', '2', 'apoptotic epicells (stats)', '2', 'dead epicells (stats)', '2', 'epicell (stats)', '2', 'expressing epicells (stats)', '2', 'incubating epicells (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotDataOverTime3Display.SeriesMarkerStyle = ['N (stats)', '0', 'Time (stats)', '0', 'apoptotic epicells (stats)', '0', 'dead epicells (stats)', '0', 'epicell (stats)', '0', 'expressing epicells (stats)', '0', 'incubating epicells (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime3Display.ShowQuartiles = 0
plotDataOverTime3Display.ShowRanges = 0

# ----------------------------------------------------------------
# setup the visualization in view 'epicellview'
# ----------------------------------------------------------------

# show data from legacyVTKReader1
legacyVTKReader1Display = Show(legacyVTKReader1, epicellview)

# get color transfer function/color map for 'epicell'
epicellLUT = GetColorTransferFunction('epicell')
epicellLUT.InterpretValuesAsCategories = 1
epicellLUT.AnnotationsInitialized = 1
epicellLUT.NanOpacity = 0.0
epicellLUT.ScalarRangeInitialized = 1.0
epicellLUT.Annotations = ['1', 'incubating', '2', 'expressing', '3', 'apoptotic', '4', 'dead']
epicellLUT.ActiveAnnotatedValues = ['1', '2', '3', '4', '1', '2', '4']
epicellLUT.IndexedColors = [0.6666666666666666, 0.6666666666666666, 0.4980392156862745, 0.6666666666666666, 0.3333333333333333, 0.0, 0.3333333333333333, 0.0, 0.0, 0.0, 0.0, 0.0]
epicellLUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0]

# get opacity transfer function/opacity map for 'epicell'
epicellPWF = GetOpacityTransferFunction('epicell')
epicellPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
legacyVTKReader1Display.Representation = 'Surface'
legacyVTKReader1Display.ColorArrayName = ['CELLS', 'epicell']
legacyVTKReader1Display.LookupTable = epicellLUT
legacyVTKReader1Display.OSPRayScaleFunction = 'PiecewiseFunction'
legacyVTKReader1Display.SelectOrientationVectors = 'None'
legacyVTKReader1Display.ScaleFactor = 5.0
legacyVTKReader1Display.SelectScaleArray = 'epicell'
legacyVTKReader1Display.GlyphType = 'Arrow'
legacyVTKReader1Display.GlyphTableIndexArray = 'epicell'
legacyVTKReader1Display.GaussianRadius = 0.25
legacyVTKReader1Display.SetScaleArray = ['POINTS', '']
legacyVTKReader1Display.ScaleTransferFunction = 'PiecewiseFunction'
legacyVTKReader1Display.OpacityArray = ['POINTS', '']
legacyVTKReader1Display.OpacityTransferFunction = 'PiecewiseFunction'
legacyVTKReader1Display.DataAxesGrid = 'GridAxesRepresentation'
legacyVTKReader1Display.PolarAxes = 'PolarAxesRepresentation'
legacyVTKReader1Display.ScalarOpacityUnitDistance = 5.210528284270439
legacyVTKReader1Display.ScalarOpacityFunction = epicellPWF
legacyVTKReader1Display.IsosurfaceValues = [0.5]

# show data from calculator2
calculator2Display = Show(calculator2, epicellview)

# get color transfer function/color map for 'virus'
virusLUT = GetColorTransferFunction('virus')
virusLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'virus'
virusPWF = GetOpacityTransferFunction('virus')
virusPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
calculator2Display.Representation = 'Points'
calculator2Display.ColorArrayName = ['CELLS', 'virus']
calculator2Display.LookupTable = virusLUT
calculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator2Display.SelectOrientationVectors = 'None'
calculator2Display.ScaleFactor = 25.0
calculator2Display.SelectScaleArray = 'virus norm'
calculator2Display.GlyphType = 'Arrow'
calculator2Display.GlyphTableIndexArray = 'virus norm'
calculator2Display.GaussianRadius = 1.25
calculator2Display.SetScaleArray = ['POINTS', '']
calculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator2Display.OpacityArray = ['POINTS', '']
calculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator2Display.DataAxesGrid = 'GridAxesRepresentation'
calculator2Display.PolarAxes = 'PolarAxesRepresentation'
calculator2Display.ScalarOpacityUnitDistance = 26.0526414213522
calculator2Display.ScalarOpacityFunction = virusPWF
calculator2Display.IsosurfaceValues = [0.047058823529411764]

# setup the color legend parameters for each legend in this view

# get color legend/bar for epicellLUT in view epicellview
epicellLUTColorBar = GetScalarBar(epicellLUT, epicellview)
epicellLUTColorBar.WindowLocation = 'AnyLocation'
epicellLUTColorBar.Position = [0.1302793274778228, 0.04683229813664577]
epicellLUTColorBar.Title = 'epicell'
epicellLUTColorBar.ComponentTitle = ''
epicellLUTColorBar.ScalarBarLength = 0.3300000000000004

# set color bar visibility
epicellLUTColorBar.Visibility = 1

# get color legend/bar for virusLUT in view epicellview
virusLUTColorBar = GetScalarBar(virusLUT, epicellview)
virusLUTColorBar.Title = 'virus'
virusLUTColorBar.ComponentTitle = ''

# set color bar visibility
virusLUTColorBar.Visibility = 1

# show color legend
legacyVTKReader1Display.SetScalarBarVisibility(epicellview, True)

# show color legend
calculator2Display.SetScalarBarVisibility(epicellview, True)

# ----------------------------------------------------------------
# setup the visualization in view 'tcellchart'
# ----------------------------------------------------------------

# show data from plotDataOverTime4
plotDataOverTime4Display = Show(plotDataOverTime4, tcellchart)

# trace defaults for the display properties.
plotDataOverTime4Display.AttributeType = 'Row Data'
plotDataOverTime4Display.UseIndexForXAxis = 0
plotDataOverTime4Display.XArrayName = 'Time'
plotDataOverTime4Display.SeriesVisibility = ['t-cell-tissue (stats)']
plotDataOverTime4Display.SeriesLabel = ['t-cell-tissue (stats)', 't-cell-tissue (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)']
plotDataOverTime4Display.SeriesColor = ['t-cell-tissue (stats)', '0', '0.666667', '0', 'N (stats)', '0.889998', '0.100008', '0.110002', 'Time (stats)', '0.220005', '0.489998', '0.719997', 'vtkValidPointMask (stats)', '0.300008', '0.689998', '0.289998']
plotDataOverTime4Display.SeriesPlotCorner = ['N (stats)', '0', 'Time (stats)', '0', 't-cell-tissue (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime4Display.SeriesLabelPrefix = ''
plotDataOverTime4Display.SeriesLineStyle = ['N (stats)', '1', 'Time (stats)', '1', 't-cell-tissue (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotDataOverTime4Display.SeriesLineThickness = ['N (stats)', '2', 'Time (stats)', '2', 't-cell-tissue (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotDataOverTime4Display.SeriesMarkerStyle = ['N (stats)', '0', 'Time (stats)', '0', 't-cell-tissue (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime4Display.ShowQuartiles = 0
plotDataOverTime4Display.ShowRanges = 0

# ----------------------------------------------------------------
# setup the visualization in view 'tcellview'
# ----------------------------------------------------------------

# show data from legacyVTKReader4
legacyVTKReader4Display = Show(legacyVTKReader4, tcellview)

# get color transfer function/color map for 'tcelltissue'
tcelltissueLUT = GetColorTransferFunction('tcelltissue')
tcelltissueLUT.InterpretValuesAsCategories = 1
tcelltissueLUT.AnnotationsInitialized = 1
tcelltissueLUT.RGBPoints = [0.0, 0.054902, 0.109804, 0.121569, 0.05, 0.07451, 0.172549, 0.180392, 0.1, 0.086275, 0.231373, 0.219608, 0.15, 0.094118, 0.278431, 0.25098, 0.2, 0.109804, 0.34902, 0.278431, 0.25, 0.113725, 0.4, 0.278431, 0.3, 0.117647, 0.45098, 0.270588, 0.35, 0.117647, 0.490196, 0.243137, 0.4, 0.113725, 0.521569, 0.203922, 0.45, 0.109804, 0.54902, 0.152941, 0.5, 0.082353, 0.588235, 0.082353, 0.55, 0.109804, 0.631373, 0.05098, 0.6, 0.211765, 0.678431, 0.082353, 0.65, 0.317647, 0.721569, 0.113725, 0.7, 0.431373, 0.760784, 0.160784, 0.75, 0.556863, 0.8, 0.239216, 0.8, 0.666667, 0.839216, 0.294118, 0.85, 0.784314, 0.878431, 0.396078, 0.9, 0.886275, 0.921569, 0.533333, 0.95, 0.960784, 0.94902, 0.670588, 1.0, 1.0, 0.984314, 0.901961]
tcelltissueLUT.ColorSpace = 'Lab'
tcelltissueLUT.NanColor = [0.25, 0.0, 0.0]
tcelltissueLUT.NanOpacity = 0.0
tcelltissueLUT.ScalarRangeInitialized = 1.0
tcelltissueLUT.Annotations = ['1', '1', '2', '2', '3', '3', '4', '4', '5', '5']
tcelltissueLUT.ActiveAnnotatedValues = ['1', '2', '4']
tcelltissueLUT.IndexedColors = [0.0, 0.9921568627450981, 0.047058823529411764, 0.43137254901960786, 0.9921568627450981, 0.0, 0.8, 1.0, 0.0, 0.9058823529411765, 1.0, 0.0392156862745098, 1.0, 0.9686274509803922, 0.0392156862745098]
tcelltissueLUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0]

# get opacity transfer function/opacity map for 'tcelltissue'
tcelltissuePWF = GetOpacityTransferFunction('tcelltissue')
tcelltissuePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
legacyVTKReader4Display.Representation = 'Surface'
legacyVTKReader4Display.ColorArrayName = ['CELLS', 't-cell-tissue']
legacyVTKReader4Display.LookupTable = tcelltissueLUT
legacyVTKReader4Display.OSPRayScaleFunction = 'PiecewiseFunction'
legacyVTKReader4Display.SelectOrientationVectors = 'None'
legacyVTKReader4Display.ScaleFactor = 5.0
legacyVTKReader4Display.SelectScaleArray = 't-cell-tissue'
legacyVTKReader4Display.GlyphType = 'Arrow'
legacyVTKReader4Display.GlyphTableIndexArray = 't-cell-tissue'
legacyVTKReader4Display.GaussianRadius = 0.25
legacyVTKReader4Display.SetScaleArray = ['POINTS', '']
legacyVTKReader4Display.ScaleTransferFunction = 'PiecewiseFunction'
legacyVTKReader4Display.OpacityArray = ['POINTS', '']
legacyVTKReader4Display.OpacityTransferFunction = 'PiecewiseFunction'
legacyVTKReader4Display.DataAxesGrid = 'GridAxesRepresentation'
legacyVTKReader4Display.PolarAxes = 'PolarAxesRepresentation'
legacyVTKReader4Display.ScalarOpacityUnitDistance = 5.210528284270439
legacyVTKReader4Display.ScalarOpacityFunction = tcelltissuePWF
legacyVTKReader4Display.IsosurfaceValues = [0.5]

# show data from calculator1
calculator1Display = Show(calculator1, tcellview)

# get color transfer function/color map for 'chemokine'
chemokineLUT = GetColorTransferFunction('chemokine')
chemokineLUT.AutomaticRescaleRangeMode = 'Never'
chemokineLUT.RGBPoints = [0.0001, 0.0, 0.0, 0.0, 1.2, 1.0, 1.0, 1.0]
chemokineLUT.ColorSpace = 'RGB'
chemokineLUT.NanColor = [1.0, 0.0, 0.0]
chemokineLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'chemokine'
chemokinePWF = GetOpacityTransferFunction('chemokine')
chemokinePWF.Points = [0.0001, 0.0, 0.5, 0.0, 1.2, 1.0, 0.5, 0.0]
chemokinePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['CELLS', 'chemokine']
calculator1Display.LookupTable = chemokineLUT
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 25.0
calculator1Display.SelectScaleArray = 'chemokine'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'chemokine'
calculator1Display.GaussianRadius = 1.25
calculator1Display.SetScaleArray = ['POINTS', '']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', '']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.ScalarOpacityUnitDistance = 26.0526414213522
calculator1Display.ScalarOpacityFunction = chemokinePWF
calculator1Display.IsosurfaceValues = [0.5]

# setup the color legend parameters for each legend in this view

# get color legend/bar for tcelltissueLUT in view tcellview
tcelltissueLUTColorBar = GetScalarBar(tcelltissueLUT, tcellview)
tcelltissueLUTColorBar.WindowLocation = 'UpperRightCorner'
tcelltissueLUTColorBar.Position = [0.03716608594657392, 0.2766084788029926]
tcelltissueLUTColorBar.Title = 't-cell-tissue'
tcelltissueLUTColorBar.ComponentTitle = ''
tcelltissueLUTColorBar.ScalarBarLength = 0.3299999999999998

# set color bar visibility
tcelltissueLUTColorBar.Visibility = 1

# get color legend/bar for chemokineLUT in view tcellview
chemokineLUTColorBar = GetScalarBar(chemokineLUT, tcellview)
chemokineLUTColorBar.Title = 'chemokine'
chemokineLUTColorBar.ComponentTitle = ''

# set color bar visibility
chemokineLUTColorBar.Visibility = 1

# show color legend
legacyVTKReader4Display.SetScalarBarVisibility(tcellview, True)

# show color legend
calculator1Display.SetScalarBarVisibility(tcellview, True)

# ----------------------------------------------------------------
# setup the visualization in view 'viruschart'
# ----------------------------------------------------------------

# show data from plotDataOverTime2
plotDataOverTime2Display = Show(plotDataOverTime2, viruschart)

# trace defaults for the display properties.
plotDataOverTime2Display.AttributeType = 'Row Data'
plotDataOverTime2Display.UseIndexForXAxis = 0
plotDataOverTime2Display.XArrayName = 'Time'
plotDataOverTime2Display.SeriesVisibility = ['virus (stats)', 'virus norm (stats)']
plotDataOverTime2Display.SeriesLabel = ['virus (stats)', 'virus (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)', 'virus norm (stats)', 'virus norm (stats)']
plotDataOverTime2Display.SeriesColor = ['virus (stats)', '0.889998', '0.100008', '0.110002', 'N (stats)', '0.220005', '0.489998', '0.719997', 'Time (stats)', '0.300008', '0.689998', '0.289998', 'vtkValidPointMask (stats)', '0.6', '0.310002', '0.639994', 'virus norm (stats)', '1', '0', '0']
plotDataOverTime2Display.SeriesPlotCorner = ['N (stats)', '0', 'Time (stats)', '0', 'virus (stats)', '0', 'virus norm (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime2Display.SeriesLabelPrefix = ''
plotDataOverTime2Display.SeriesLineStyle = ['N (stats)', '1', 'Time (stats)', '1', 'virus (stats)', '1', 'virus norm (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotDataOverTime2Display.SeriesLineThickness = ['N (stats)', '2', 'Time (stats)', '2', 'virus (stats)', '2', 'virus norm (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotDataOverTime2Display.SeriesMarkerStyle = ['N (stats)', '0', 'Time (stats)', '0', 'virus (stats)', '0', 'virus norm (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime2Display.ShowQuartiles = 0
plotDataOverTime2Display.ShowRanges = 0

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------

##### INSERT BEGIN
state_fname = options.data + '-state.pvsm'
servermanager.SaveState(state_fname)

print('Created a state file', state_fname, 'for paraview')
print('Run with\n', 'paraview', state_fname)
##### INSERT END
