#!/usr/bin/env pvpython

import glob
import os
import sys
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
    render_view.CameraPosition = [xdim / 2, ydim / 2, 1480]
    render_view.CameraFocalPoint = [xdim / 2, ydim / 2, 2.5]
    return render_view


def get_sample_fnames(prefix, data_dir):
    return sorted(glob.glob(data_dir + '/' + prefix + '*'), key=lambda s: int(s.split('_')[2].split('.')[0]))


def display_data(prefix, data_dir, render_view, label, color_func):
    reader = pvs.LegacyVTKReader(FileNames=get_sample_fnames(prefix, data_dir))
    display = pvs.Show(reader, render_view)
    display.Representation = 'Surface'
    display.ColorArrayName = ['CELLS', label]
    display.LookupTable = color_func
    display.ScaleFactor = 50.0
    display.SelectScaleArray = label
    display.GlyphTableIndexArray = label
    display.UseSeparateColorMap = True
    color_bar = pvs.GetScalarBar(color_func, render_view)
    color_bar.Title = label
    color_bar.ComponentTitle = ''
    color_bar.ScalarBarLength = 0.12
    display.SetScalarBarVisibility(render_view, True)

    
def main():   
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-d', '--data', required=True, help='Specify the directory containing the time series data')
    argparser.add_argument('-o', '--output', help='The name for the state file to be loaded by paraview, without an extension')
    options = argparser.parse_args()

    xdim, ydim = get_dims(options.data)
    
    virions_view = create_render_view('virions', xdim, ydim)
    epicells_view = create_render_view('epicells', xdim, ydim)
    tcells_view = create_render_view('chemokines and tcells', xdim, ydim)

    pvs.AddCameraLink(virions_view, epicells_view, 'link1')
    pvs.AddCameraLink(epicells_view, tcells_view, 'link2')

    virions_color_func = pvs.GetColorTransferFunction('virus')
    virions_color_func.RGBPoints = [0.0255, 0.231, 0.298, 0.753, 2.55, 0.865, 0.865, 0.865, 255.0, 0.706, 0.0157, 0.149]
    virions_color_func.UseLogScale = 1
    display_data('sample_virus_', options.data, virions_view, 'virus', virions_color_func)
    
    epicells_color_func = pvs.GetColorTransferFunction('epicell')
    epicells_color_func.InterpretValuesAsCategories = 1
    epicells_color_func.Annotations = ['1', 'healthy', '2', 'incubating', '3', 'expressing', '4', 'apoptotic', '5', 'dead']
    epicells_color_func.IndexedColors = [0.8, 0.8, 0.8, 1.0, 0.333, 1.0, 0.667, 0.0, 0.0, 1.0, 1.0, 0.0, 0.333, 0.333, 0.0]
    display_data('sample_epicell_', options.data, epicells_view, 'epicell', epicells_color_func)

    chemokine_color_func = pvs.GetColorTransferFunction('chemokine')
    chemokine_color_func.RGBPoints = [0.0249, 1.0, 1.0, 1.0, 249.0, 0.0, 0.0, 0.0]
    chemokine_color_func.UseLogScale = 1
    chemokine_color_func.NanOpacity = 0.0
    display_data('sample_chemokine_', options.data, tcells_view, 'chemokine', chemokine_color_func)

    tcells_color_func = pvs.GetColorTransferFunction('t-cell-tissue')
    tcells_color_func.InterpretValuesAsCategories = 1
    tcells_color_func.NanOpacity = 0.0
    tcells_color_func.Annotations = ['255', '']
    tcells_color_func.IndexedColors = [0.0, 1.0, 0.0]
    display_data('sample_tcelltissue_', options.data, tcells_view, 't-cell-tissue', tcells_color_func)

    layout = pvs.CreateLayout(name='Main Layout')
    layout.SplitHorizontal(0, 0.33333)
    layout.SplitHorizontal(2, 0.500000)
    layout.AssignView(1, virions_view)
    layout.AssignView(5, epicells_view)
    layout.AssignView(6, tcells_view)

    if not options.output:
        state_fname = options.data + '-state.pvsm'
    else:
        state_fname = options.output + '.pvsm'

    pvs.servermanager.SaveState(state_fname)
    print('Created a state file', state_fname, 'for paraview')
    print('Run with\n', 'paraview', state_fname)


    
if __name__ == "__main__":
    main()


    







