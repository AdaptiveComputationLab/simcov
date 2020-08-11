# state file generated using paraview version 5.7.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

import glob
import os
import sys
import argparse

#### import the simple module from the paraview
from paraview.simple import *

argparser = argparse.ArgumentParser()
argparser.add_argument("--data", required=True, help='Specify the output directory containing the time series data')
options = argparser.parse_args()

def get_sample_fnames(prefix):
    global options
    num_fnames = len(glob.glob(options.data + '/' + prefix + '*'))
    return [options.data + '/' + prefix + str(i) + '.vtk' for i in range(num_fnames)]

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Quartile Chart View'
quartileChartView1 = CreateView('QuartileChartView')
quartileChartView1.ViewSize = [555, 226]
quartileChartView1.LegendPosition = [400, 501]
quartileChartView1.LeftAxisRangeMaximum = 2500.0
quartileChartView1.BottomAxisRangeMaximum = 200.0
quartileChartView1.RightAxisRangeMaximum = 6.66
quartileChartView1.TopAxisRangeMaximum = 6.66

# Create a new 'Quartile Chart View'
quartileChartView2 = CreateView('QuartileChartView')
quartileChartView2.ViewSize = [555, 226]
quartileChartView2.LegendPosition = [425, 501]
quartileChartView2.LeftAxisRangeMaximum = 250000.0
quartileChartView2.BottomAxisRangeMaximum = 200.0
quartileChartView2.RightAxisRangeMaximum = 6.66
quartileChartView2.TopAxisRangeMaximum = 6.66

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1120, 901]
renderView1.InteractionMode = 'Selection'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [100.0, 100.0, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [100.0, 100.0, 702.5932558000225]
renderView1.CameraFocalPoint = [100.0, 100.0, 0.5]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 42.875511494162964
renderView1.Background = [0.32, 0.34, 0.43]

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [1120, 644]
renderView2.InteractionMode = 'Selection'
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [100.0, 100.0, 0.5]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [100.0, 100.0, 702.5932558000225]
renderView2.CameraFocalPoint = [100.0, 100.0, 0.5]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 42.875511494162964
renderView2.Background = [0.32, 0.34, 0.43]

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.SplitVertical(2, 0.725764)
layout1.AssignView(5, renderView2)
layout1.SplitHorizontal(6, 0.500000)
layout1.AssignView(13, quartileChartView1)
layout1.AssignView(14, quartileChartView2)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
chemokine_reader = LegacyVTKReader(FileNames=get_sample_fnames('sample_chemokine_'))

# create a new 'Calculator'
calculator3 = Calculator(Input=chemokine_reader)
calculator3.AttributeType = 'Cell Data'
calculator3.ResultArrayName = 'chemokine norm'
calculator3.Function = 'chemokine/255'

# create a new 'Legacy VTK Reader'
tcell_reader = LegacyVTKReader(FileNames=get_sample_fnames('sample_tcelltissue_'))

# create a new 'Legacy VTK Reader'
virus_reader = LegacyVTKReader(FileNames=get_sample_fnames('sample_virus_'))

# create a new 'Calculator'
calculator2 = Calculator(Input=virus_reader)
calculator2.AttributeType = 'Cell Data'
calculator2.ResultArrayName = 'virus norm'
calculator2.Function = 'virus/255'

# create a new 'Legacy VTK Reader'
icytokine_reader = LegacyVTKReader(FileNames=get_sample_fnames('sample_icytokine_'))

# create a new 'Calculator'
calculator1 = Calculator(Input=icytokine_reader)
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'icytokine norm'
calculator1.Function = 'icytokine/255'

# create a new 'Integrate Variables'
integrateVariables2 = IntegrateVariables(Input=tcell_reader)

# create a new 'Plot Data Over Time'
plotDataOverTime1 = PlotDataOverTime(Input=integrateVariables2)
plotDataOverTime1.FieldAssociation = 'Cells'

# create a new 'Legacy VTK Reader'
epicell_reader = LegacyVTKReader(FileNames=get_sample_fnames('sample_epicell_'))

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=epicell_reader)

# create a new 'Plot Data Over Time'
plotDataOverTime2 = PlotDataOverTime(Input=integrateVariables1)
plotDataOverTime2.FieldAssociation = 'Cells'

# ----------------------------------------------------------------
# setup the visualization in view 'quartileChartView1'
# ----------------------------------------------------------------

# show data from plotDataOverTime1
plotDataOverTime1Display = Show(plotDataOverTime1, quartileChartView1)

# trace defaults for the display properties.
plotDataOverTime1Display.AttributeType = 'Row Data'
plotDataOverTime1Display.UseIndexForXAxis = 0
plotDataOverTime1Display.XArrayName = 'Time'
plotDataOverTime1Display.SeriesVisibility = ['t-cell-tissue (stats)']
plotDataOverTime1Display.SeriesLabel = ['t-cell-tissue (stats)', 't-cell-tissue (stats)', 'Volume (stats)', 'Volume (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)']
plotDataOverTime1Display.SeriesColor = ['t-cell-tissue (stats)', '0', '0.666667', '0', 'Volume (stats)', '0.889998', '0.100008', '0.110002', 'N (stats)', '0.220005', '0.489998', '0.719997', 'Time (stats)', '0.300008', '0.689998', '0.289998', 'vtkValidPointMask (stats)', '0.6', '0.310002', '0.639994']
plotDataOverTime1Display.SeriesPlotCorner = ['N (stats)', '0', 'Time (stats)', '0', 'Volume (stats)', '0', 't-cell-tissue (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime1Display.SeriesLabelPrefix = ''
plotDataOverTime1Display.SeriesLineStyle = ['N (stats)', '1', 'Time (stats)', '1', 'Volume (stats)', '1', 't-cell-tissue (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotDataOverTime1Display.SeriesLineThickness = ['N (stats)', '2', 'Time (stats)', '2', 'Volume (stats)', '2', 't-cell-tissue (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotDataOverTime1Display.SeriesMarkerStyle = ['N (stats)', '0', 'Time (stats)', '0', 'Volume (stats)', '0', 't-cell-tissue (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime1Display.ShowQuartiles = 0
plotDataOverTime1Display.ShowRanges = 0

# ----------------------------------------------------------------
# setup the visualization in view 'quartileChartView2'
# ----------------------------------------------------------------

# show data from plotDataOverTime2
plotDataOverTime2Display = Show(plotDataOverTime2, quartileChartView2)

# trace defaults for the display properties.
plotDataOverTime2Display.AttributeType = 'Row Data'
plotDataOverTime2Display.UseIndexForXAxis = 0
plotDataOverTime2Display.XArrayName = 'Time'
plotDataOverTime2Display.SeriesVisibility = ['epicell (stats)']
plotDataOverTime2Display.SeriesLabel = ['epicell (stats)', 'epicell (stats)', 'Volume (stats)', 'Volume (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)']
plotDataOverTime2Display.SeriesColor = ['epicell (stats)', '0', '0', '0', 'Volume (stats)', '0.889998', '0.100008', '0.110002', 'N (stats)', '0.220005', '0.489998', '0.719997', 'Time (stats)', '0.300008', '0.689998', '0.289998', 'vtkValidPointMask (stats)', '0.6', '0.310002', '0.639994']
plotDataOverTime2Display.SeriesPlotCorner = ['N (stats)', '0', 'Time (stats)', '0', 'Volume (stats)', '0', 'epicell (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime2Display.SeriesLabelPrefix = ''
plotDataOverTime2Display.SeriesLineStyle = ['N (stats)', '1', 'Time (stats)', '1', 'Volume (stats)', '1', 'epicell (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotDataOverTime2Display.SeriesLineThickness = ['N (stats)', '2', 'Time (stats)', '2', 'Volume (stats)', '2', 'epicell (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotDataOverTime2Display.SeriesMarkerStyle = ['N (stats)', '0', 'Time (stats)', '0', 'Volume (stats)', '0', 'epicell (stats)', '0', 'vtkValidPointMask (stats)', '0']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from epicell_reader
epicell_readerDisplay = Show(epicell_reader, renderView1)

# get color transfer function/color map for 'epicell'
epicellLUT = GetColorTransferFunction('epicell')
epicellLUT.InterpretValuesAsCategories = 1
epicellLUT.AnnotationsInitialized = 1
epicellLUT.NanOpacity = 0.0
epicellLUT.ScalarRangeInitialized = 1.0
epicellLUT.Annotations = ['1', 'iincubating', '2', 'expressing', '3', 'apoptotic', '4', 'dead']
epicellLUT.ActiveAnnotatedValues = ['1', '2', '3', '4', '1', '2', '4']
epicellLUT.IndexedColors = [0.6666666666666666, 0.6666666666666666, 0.4980392156862745, 0.6666666666666666, 0.3333333333333333, 0.0, 0.3333333333333333, 0.0, 0.0, 0.0, 0.0, 0.0]
epicellLUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0]

# get opacity transfer function/opacity map for 'epicell'
epicellPWF = GetOpacityTransferFunction('epicell')
epicellPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
epicell_readerDisplay.Representation = 'Surface With Edges'
epicell_readerDisplay.ColorArrayName = ['CELLS', 'epicell']
epicell_readerDisplay.LookupTable = epicellLUT
epicell_readerDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
epicell_readerDisplay.SelectOrientationVectors = 'None'
epicell_readerDisplay.ScaleFactor = 5.0
epicell_readerDisplay.SelectScaleArray = 'epicell'
epicell_readerDisplay.GlyphType = 'Arrow'
epicell_readerDisplay.GlyphTableIndexArray = 'epicell'
epicell_readerDisplay.GaussianRadius = 0.25
epicell_readerDisplay.SetScaleArray = ['POINTS', '']
epicell_readerDisplay.ScaleTransferFunction = 'PiecewiseFunction'
epicell_readerDisplay.OpacityArray = ['POINTS', '']
epicell_readerDisplay.OpacityTransferFunction = 'PiecewiseFunction'
epicell_readerDisplay.DataAxesGrid = 'GridAxesRepresentation'
epicell_readerDisplay.PolarAxes = 'PolarAxesRepresentation'
epicell_readerDisplay.ScalarOpacityUnitDistance = 5.210528284270439
epicell_readerDisplay.ScalarOpacityFunction = epicellPWF
epicell_readerDisplay.IsosurfaceValues = [0.5]

# setup the color legend parameters for each legend in this view

# get color legend/bar for epicellLUT in view renderView1
epicellLUTColorBar = GetScalarBar(epicellLUT, renderView1)
epicellLUTColorBar.Orientation = 'Horizontal'
epicellLUTColorBar.WindowLocation = 'AnyLocation'
epicellLUTColorBar.Position = [0.015535714285714236, 0.038855098389982246]
epicellLUTColorBar.Title = 'epicell'
epicellLUTColorBar.ComponentTitle = ''
epicellLUTColorBar.ScalarBarLength = 0.3300000000000003

# set color bar visibility
epicellLUTColorBar.Visibility = 1

# show color legend
epicell_readerDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from tcell_reader
tcell_readerDisplay = Show(tcell_reader, renderView2)

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
tcell_readerDisplay.Representation = 'Surface With Edges'
tcell_readerDisplay.ColorArrayName = ['CELLS', 't-cell-tissue']
tcell_readerDisplay.LookupTable = tcelltissueLUT
tcell_readerDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
tcell_readerDisplay.SelectOrientationVectors = 'None'
tcell_readerDisplay.ScaleFactor = 5.0
tcell_readerDisplay.SelectScaleArray = 't-cell-tissue'
tcell_readerDisplay.GlyphType = 'Arrow'
tcell_readerDisplay.GlyphTableIndexArray = 't-cell-tissue'
tcell_readerDisplay.GaussianRadius = 0.25
tcell_readerDisplay.SetScaleArray = ['POINTS', '']
tcell_readerDisplay.ScaleTransferFunction = 'PiecewiseFunction'
tcell_readerDisplay.OpacityArray = ['POINTS', '']
tcell_readerDisplay.OpacityTransferFunction = 'PiecewiseFunction'
tcell_readerDisplay.DataAxesGrid = 'GridAxesRepresentation'
tcell_readerDisplay.PolarAxes = 'PolarAxesRepresentation'
tcell_readerDisplay.ScalarOpacityUnitDistance = 5.210528284270439
tcell_readerDisplay.ScalarOpacityFunction = tcelltissuePWF
tcell_readerDisplay.IsosurfaceValues = [0.5]

# show data from calculator3
calculator3Display = Show(calculator3, renderView2)

# get color transfer function/color map for 'chemokinenorm'
chemokinenormLUT = GetColorTransferFunction('chemokinenorm')
chemokinenormLUT.RGBPoints = [0.0001, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
chemokinenormLUT.UseLogScale = 1
chemokinenormLUT.ColorSpace = 'RGB'
chemokinenormLUT.NanColor = [1.0, 0.0, 0.0]
chemokinenormLUT.NanOpacity = 0.0
chemokinenormLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'chemokinenorm'
chemokinenormPWF = GetOpacityTransferFunction('chemokinenorm')
chemokinenormPWF.Points = [0.0001, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
chemokinenormPWF.UseLogScale = 1
chemokinenormPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
calculator3Display.Representation = 'Surface With Edges'
calculator3Display.ColorArrayName = ['CELLS', 'chemokine norm']
calculator3Display.LookupTable = chemokinenormLUT
calculator3Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator3Display.SelectOrientationVectors = 'None'
calculator3Display.ScaleFactor = 25.0
calculator3Display.SelectScaleArray = 'chemokine'
calculator3Display.GlyphType = 'Arrow'
calculator3Display.GlyphTableIndexArray = 'chemokine'
calculator3Display.GaussianRadius = 1.25
calculator3Display.SetScaleArray = ['POINTS', '']
calculator3Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator3Display.OpacityArray = ['POINTS', '']
calculator3Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator3Display.DataAxesGrid = 'GridAxesRepresentation'
calculator3Display.PolarAxes = 'PolarAxesRepresentation'
calculator3Display.ScalarOpacityUnitDistance = 26.0526414213522
calculator3Display.ScalarOpacityFunction = chemokinenormPWF
calculator3Display.IsosurfaceValues = [0.5]

# setup the color legend parameters for each legend in this view

# get color legend/bar for tcelltissueLUT in view renderView2
tcelltissueLUTColorBar = GetScalarBar(tcelltissueLUT, renderView2)
tcelltissueLUTColorBar.WindowLocation = 'UpperRightCorner'
tcelltissueLUTColorBar.Position = [0.03716608594657392, 0.2766084788029926]
tcelltissueLUTColorBar.Title = 't-cell-tissue'
tcelltissueLUTColorBar.ComponentTitle = ''
tcelltissueLUTColorBar.ScalarBarLength = 0.3299999999999998

# set color bar visibility
tcelltissueLUTColorBar.Visibility = 1

# get color legend/bar for chemokinenormLUT in view renderView2
chemokinenormLUTColorBar = GetScalarBar(chemokinenormLUT, renderView2)
chemokinenormLUTColorBar.WindowLocation = 'AnyLocation'
chemokinenormLUTColorBar.Position = [0.020905923344947744, 0.6456857855361595]
chemokinenormLUTColorBar.Title = 'chemokine norm'
chemokinenormLUTColorBar.ComponentTitle = ''
chemokinenormLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
chemokinenormLUTColorBar.Visibility = 1

# show color legend
tcell_readerDisplay.SetScalarBarVisibility(renderView2, True)

# show color legend
calculator3Display.SetScalarBarVisibility(renderView2, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------

state_fname = os.path.splitext(sys.argv[0])[0] + '.pvsm'
servermanager.SaveState(state_fname)

print('Created a state file', state_fname, 'for paraview')
print('Run with\n', 'paraview', state_fname)
