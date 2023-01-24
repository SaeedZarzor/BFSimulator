# trace generated using paraview version 5.10.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
from pathlib import Path
import os
import re
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

directory_folder = "Folder_Output"
parent_dir_folder = os.getcwd()
path_folder = os.path.join(parent_dir_folder, directory_folder)
Files = Path(path_folder)
# Iterate directory
res=[]
a=0.0
for f in Files.iterdir():
    for file in Files.iterdir():
        if re.match("Output_3_\d+\.vtk", file.name):
            Output, degree, number=file.stem.split('_')
            if(a == float(number)):
                res.append(str(file))
                a+=1
# create a new 'Legacy VTK Reader'
output_3_ = LegacyVTKReader(registrationName='Output_3_*', FileNames= res)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
output_2_Display = Show(output_3_, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'cell_desnity'
cell_desnityLUT = GetColorTransferFunction('cell_desnity')

# get opacity transfer function/opacity map for 'cell_desnity'
cell_desnityPWF = GetOpacityTransferFunction('cell_desnity')

# trace defaults for the display properties.
output_2_Display.Representation = 'Surface'
output_2_Display.ColorArrayName = ['POINTS', 'cell_desnity']
output_2_Display.LookupTable = cell_desnityLUT
output_2_Display.SelectTCoordArray = 'None'
output_2_Display.SelectNormalArray = 'None'
output_2_Display.SelectTangentArray = 'None'
output_2_Display.OSPRayScaleArray = 'cell_desnity'
output_2_Display.OSPRayScaleFunction = 'PiecewiseFunction'
output_2_Display.SelectOrientationVectors = 'displacement'
output_2_Display.ScaleFactor = 0.4
output_2_Display.SelectScaleArray = 'cell_desnity'
output_2_Display.GlyphType = 'Arrow'
output_2_Display.GlyphTableIndexArray = 'cell_desnity'
output_2_Display.GaussianRadius = 0.02
output_2_Display.SetScaleArray = ['POINTS', 'cell_desnity']
output_2_Display.ScaleTransferFunction = 'PiecewiseFunction'
output_2_Display.OpacityArray = ['POINTS', 'cell_desnity']
output_2_Display.OpacityTransferFunction = 'PiecewiseFunction'
output_2_Display.DataAxesGrid = 'GridAxesRepresentation'
output_2_Display.PolarAxes = 'PolarAxesRepresentation'
output_2_Display.ScalarOpacityFunction = cell_desnityPWF
output_2_Display.ScalarOpacityUnitDistance = 0.43868963652317955
output_2_Display.OpacityArrayName = ['POINTS', 'cell_desnity']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
output_2_Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
output_2_Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

# show color bar/color legend
output_2_Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(output_2_Display, ('POINTS', 'Stiffness', 'X'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(cell_desnityLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
output_2_Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
output_2_Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Stiffness'
stiffnessLUT = GetColorTransferFunction('Stiffness')

# get opacity transfer function/opacity map for 'Stiffness'
stiffnessPWF = GetOpacityTransferFunction('Stiffness')

animationScene1.GoToLast()

# hide color bar/color legend
output_2_Display.SetScalarBarVisibility(renderView1, False)

# show color bar/color legend
output_2_Display.SetScalarBarVisibility(renderView1, True)

# get color legend/bar for stiffnessLUT in view renderView1
stiffnessLUTColorBar = GetScalarBar(stiffnessLUT, renderView1)

# Rescale transfer function
stiffnessLUT.RescaleTransferFunction(0.0, 2.569)

# Rescale transfer function
stiffnessPWF.RescaleTransferFunction(0.0, 2.569)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(2000, 1354)

# current camera placement for renderView1
renderView1.CameraPosition = [2.0038344494845637, 7.835589358648983, 9.145605621453992]
renderView1.CameraFocalPoint = [9.46977582464255e-17, 0.9982166949999993, -4.90044343373111e-17]
renderView1.CameraViewUp = [-0.14170690174329095, 0.8074741458843637, -0.5726295990661283]
renderView1.CameraParallelScale = 3.0005949060439203

# save animation
SaveAnimation(str(path_folder)+'/Stiffness.avi', renderView1, ImageResolution=[2000, 1352],
    OverrideColorPalette='WhiteBackground',
    FrameRate=10,
    FrameWindow=[0, len(res)])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(2000, 1354)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [2.0038344494845637, 7.835589358648983, 9.145605621453992]
renderView1.CameraFocalPoint = [9.46977582464255e-17, 0.9982166949999993, -4.90044343373111e-17]
renderView1.CameraViewUp = [-0.14170690174329095, 0.8074741458843637, -0.5726295990661283]
renderView1.CameraParallelScale = 3.0005949060439203

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
