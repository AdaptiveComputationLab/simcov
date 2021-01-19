#!/usr/bin/env pvpython

import glob
import os
import sys
import math
import argparse
import paraview.simple as pvs

def get_dims(data_dir):
    xdim = 0
    ydim = 0
    xspace = 0
    yspace = 0
    with open(data_dir + '/sample_epicell_0.vtk') as f:
        for line in f.readlines():
            if line.startswith('DIMENSIONS'):
                xdim, ydim, _ = map(int, line.split()[1:])
                print('Dimensions of simulation are:', xdim - 1, ydim - 1)
            if line.startswith('SPACING'):
                xspace, yspace, _ = map(int, line.split()[1:])
                print('Spacings of simulation are:', xspace, yspace)
                break
        else:
            print('Cannot find dimensions in sample files', file=sys.stderr)
            os.abort()
    return (xdim - 1) * xspace, (ydim - 1) * yspace


def create_render_view(label, xdim, ydim):
    render_view = pvs.CreateView('RenderView')
    pvs.RenameView(label, render_view)
    render_view.OrientationAxesVisibility = 0
    render_view.CameraPosition = [xdim / 2, ydim / 2, 10000]
    render_view.CameraFocalPoint = [xdim / 2, ydim / 2, 2.5]
    render_view.AxesGrid.Visibility = 1
    #print(dir(render_view))
    return render_view


def create_chart_view(ylabel, log_scale=False, left_max_range=1000):
    chart_view = pvs.CreateView('XYChartView')
    pvs.RenameView(ylabel, chart_view)
    chart_view.LeftAxisTitle = ylabel
    chart_view.LeftAxisTitleFontSize = 12
    chart_view.LeftAxisGridColor = [0.7, 0.7, 0.7]
    chart_view.BottomAxisTitle = 'Time step'
    chart_view.BottomAxisTitleFontSize = 12
    chart_view.BottomAxisGridColor = [0.7, 0.7, 0.7]
    chart_view.BottomAxisUseCustomRange = 0
    #print(dir(chart_view))
    chart_view.LegendLocation = 'TopLeft'
    if log_scale:
        chart_view.LeftAxisLogScale = 1
        chart_view.LeftAxisUseCustomRange = 1
        chart_view.LeftAxisRangeMinimum = 1.0
        chart_view.LeftAxisRangeMaximum = left_max_range
    return chart_view

    
def get_sample_fnames(prefix, data_dir):
    return sorted(glob.glob(data_dir + '/' + prefix + '*'), key=lambda s: int(s.split('_')[2].split('.')[0]))


def display_data(prefix, data_dir, render_view, label, color_func, representation='Surface'):
    reader = pvs.LegacyVTKReader(FileNames=get_sample_fnames(prefix, data_dir), registrationName=label)
    display = pvs.Show(reader, render_view)
    display.Representation = representation
    display.ColorArrayName = ['CELLS', label]
    display.LookupTable = color_func
    display.ScaleFactor = 50.0
    display.SelectScaleArray = label
    display.GlyphTableIndexArray = label
    display.UseSeparateColorMap = True
    #print(dir(display))
    color_bar = pvs.GetScalarBar(color_func, render_view)
    color_bar.Title = label
    color_bar.ComponentTitle = ''
    color_bar.ScalarBarLength = 0.12
    display.SetScalarBarVisibility(render_view, True)


def display_chart(stats_fname, chart_view, cols, col_labels, label, colors=None):
    reader = pvs.CSVReader(FileName=[stats_fname], registrationName=label + 'Reader')
    reader.AddTabFieldDelimiter = 1
    # for a plot of data over time - don't know how to get it to work with CSVReader data - missing time series
    # plot_data = pvs.PlotDataOverTime(Input=reader, registrationName=label)
    # plot_data.FieldAssociation = 'Cells'
    # plot_display = pvs.Show(plot_data, chart_view)
    plot_display = pvs.Show(reader, chart_view)
    plot_display.CompositeDataSetIndex = [0]
    plot_display.AttributeType = 'Row Data'
    plot_display.SeriesVisibility = cols
    plot_display.SeriesLabel = col_labels
    if colors is not None:
        plot_display.SeriesColor = colors
    #print(dir(plot_display))

    
def main():   
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-d', '--data', required=True, help='Specify the directory containing the time series data')
    argparser.add_argument('-s', '--stats', help='The file name containing global statistics')
    argparser.add_argument('-o', '--output', help='The name for the state file to be loaded by paraview, without an extension')
    options = argparser.parse_args()

    xdim, ydim = get_dims(options.data)
    
    virions_view = create_render_view('virions', xdim, ydim)
    virions_chart_view = create_chart_view('Average virions')
    epicells_view = create_render_view('epicells', xdim, ydim)
    epicells_chart_view = create_chart_view('epicells', log_scale=True, left_max_range=xdim * ydim)
    tcells_view = create_render_view('chemokines and tcells', xdim, ydim)
    tcells_chart_view = create_chart_view('tcells', log_scale=True, left_max_range=xdim * ydim)
    
    pvs.AddCameraLink(virions_view, epicells_view, 'link1')
    pvs.AddCameraLink(epicells_view, tcells_view, 'link2')

    virions_color_func = pvs.GetColorTransferFunction('virus')
    virions_color_func.RGBPoints = [0.0255, 0.231, 0.298, 0.753, 2.55, 0.865, 0.865, 0.865, 255.0, 0.706, 0.0157, 0.149]
    virions_color_func.UseLogScale = 0
    display_data('sample_virus_', options.data, virions_view, 'virus', virions_color_func)
    display_chart(options.stats, virions_chart_view, ['virs'], ['virs', 'virions'], 'virions', ['virs', '1', '0', '0'])

    epicell_cols = ['incb', '1.0', '0.333', '1.0', 'expr', '0.667', '0.0', '0.0', 'apop', '1.0', '1.0', '0',
                    'dead', '0.333', '0.333', '0.0']
    epicell_cols_vals = [0.8, 0.8, 0.8]
    for i in range(len(epicell_cols)):
        if i % 4 == 0: continue
        epicell_cols_vals.append(float(epicell_cols[i]))
    
    epicells_color_func = pvs.GetColorTransferFunction('epicell')
    epicells_color_func.InterpretValuesAsCategories = 1
    epicells_color_func.Annotations = ['1', 'healthy', '2', 'incubating', '3', 'expressing', '4', 'apoptotic', '5', 'dead']
    epicells_color_func.IndexedColors = epicell_cols_vals 
    
    display_data('sample_epicell_', options.data, epicells_view, 'epicell', epicells_color_func)
    display_chart(options.stats, epicells_chart_view, ['incb', 'expr', 'apop', 'dead'],
                  ['incb', 'incubating', 'expr', 'expressing', 'apop', 'apoptotic', 'dead', 'dead'], 'epicells', epicell_cols)
    
    chemokine_color_func = pvs.GetColorTransferFunction('chemokine')
    chemokine_color_func.RGBPoints = [0.0249, 0.0, 0.0, 0.0, 249.0, 1.0, 1.0, 1.0]
    chemokine_color_func.UseLogScale = 0
    chemokine_color_func.NanOpacity = 0.0
    display_data('sample_chemokine_', options.data, tcells_view, 'chemokine', chemokine_color_func)

    tcells_color_func = pvs.GetColorTransferFunction('t-cell-tissue')
    tcells_color_func.InterpretValuesAsCategories = 1
    tcells_color_func.NanOpacity = 0.0
    tcells_color_func.Annotations = ['1', '>0', '2', '>0.12', '3', '>0.25', '4', '>0.5']
    tcells_color_func.IndexedColors = [0.0, 0.3, 0.0, 0.0, 0.5, 0.0, 0.0, 0.7, 0.0, 0.0, 1.0, 0.0]
    display_data('sample_tcelltissue_', options.data, tcells_view, 't-cell-tissue', tcells_color_func)#, representation='Points')
    display_chart(options.stats, tcells_chart_view, ['ttis', 'tvas'], ['ttis', 'tissue', 'tvas', 'vasculature'], 'tcells',
                  ['ttis', '0', '1', '0', 'tvas', '0', '0', '1'])

    layout = pvs.CreateLayout(name='Main Layout')
    layout.SplitHorizontal(0, 0.333)
    layout.SplitVertical(1, 0.7)
    layout.SplitHorizontal(2, 0.5)
    layout.SplitVertical(5, 0.7)
    layout.SplitVertical(6, 0.7)
    layout.AssignView(3, virions_view)
    layout.AssignView(4, virions_chart_view)
    layout.AssignView(11, epicells_view)
    layout.AssignView(12, epicells_chart_view)
    layout.AssignView(13, tcells_view)
    layout.AssignView(14, tcells_chart_view)

    # for view in [virions_view, epicells_view, tcells_view]:
    #      print(view.GetActiveCamera().GetViewUp())
    #      print(view.GetActiveCamera().GetPosition())
    #      print(view.GetActiveCamera().GetFocalPoint())
    # for view in [virions_view, epicells_view, tcells_view]:
    #     camera = view.GetActiveCamera()
    #     focal_point = camera.GetFocalPoint()
    #     pos = camera.GetPosition()
    #     dist = math.sqrt((pos[0] - focal_point[0]) ** 2 + (pos[1] - focal_point[1]) ** 2 + (pos[2] - focal_point[2]) ** 2)
    #     camera.SetPosition(focal_point[0], focal_point[1], focal_point[2] + dist)
    #     camera.SetViewUp(0.0, 1.0, 0.0)
    
    if not options.output:
        state_fname = options.data + '-state.pvsm'
    else:
        state_fname = options.output + '.pvsm'

    pvs.servermanager.SaveState(state_fname)
    print('Created a state file', state_fname, 'for paraview')
    print('Run with\n', 'paraview', state_fname)


    
if __name__ == "__main__":
    main()


    







