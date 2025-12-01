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

animationScene1.GoToLast()

# set scalar coloring
ColorBy(output_2_Display, ('POINTS', 'growth_factor', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(cell_desnityLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
output_2_Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
output_2_Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'growth_factor'
growth_factorLUT = GetColorTransferFunction('growth_factor')

# get opacity transfer function/opacity map for 'growth_factor'
growth_factorPWF = GetOpacityTransferFunction('growth_factor')

# set scalar coloring
ColorBy(output_2_Display, ('POINTS', 'growth_factor', 'X'))

# rescale color and/or opacity maps used to exactly fit the current data range
output_2_Display.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(growth_factorLUT, output_2_Display)

# Rescale transfer function
growth_factorLUT.RescaleTransferFunction(1.29128, 1.973)

# Rescale transfer function
growth_factorPWF.RescaleTransferFunction(1.29128, 1.973)

# Rescale transfer function
growth_factorLUT.RescaleTransferFunction(0.0, 2.0)

# Rescale transfer function
growth_factorPWF.RescaleTransferFunction(0.0, 2.0)

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=output_3_)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'cell_desnity']
clip1.Value = 2764.2785

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.06186499999999984, 1.128086695, 0.0015650000000000386]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [0.06186499999999984, 1.128086695, 0.0015650000000000386]

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['POINTS', 'cell_desnity']
clip1Display.LookupTable = cell_desnityLUT
clip1Display.SelectTCoordArray = 'None'
clip1Display.SelectNormalArray = 'None'
clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleArray = 'cell_desnity'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'displacement'
clip1Display.ScaleFactor = 0.4465410000000001
clip1Display.SelectScaleArray = 'cell_desnity'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'cell_desnity'
clip1Display.GaussianRadius = 0.02232705
clip1Display.SetScaleArray = ['POINTS', 'cell_desnity']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'cell_desnity']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = cell_desnityPWF
clip1Display.ScalarOpacityUnitDistance = 0.4695967559792564
clip1Display.OpacityArrayName = ['POINTS', 'cell_desnity']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [-519.211, 0.0, 0.5, 0.0, 5657.12, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [-519.211, 0.0, 0.5, 0.0, 5657.12, 1.0, 0.5, 0.0]

# hide data in view
Hide(output_3_, renderView1)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'growth_factor', 'X'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(cell_desnityLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
growth_factorLUT.RescaleTransferFunction(0.9675437496193561, 1.93895)

# Rescale transfer function
growth_factorPWF.RescaleTransferFunction(0.9675437496193561, 1.93895)

# Rescale transfer function
growth_factorLUT.RescaleTransferFunction(1.0, 1.93895)

# Rescale transfer function
growth_factorPWF.RescaleTransferFunction(1.0, 1.93895)

# get color legend/bar for growth_factorLUT in view renderView1
growth_factorLUTColorBar = GetScalarBar(growth_factorLUT, renderView1)

# Properties modified on growth_factorLUTColorBar
growth_factorLUTColorBar.TitleColor = [0.06465247577630275, 0.11065842679484245, 0.04505989166094453]
growth_factorLUTColorBar.LabelColor = [0.06465247577630275, 0.11065842679484245, 0.04505989166094453]

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(2000, 1354)

# current camera placement for renderView1
renderView1.CameraPosition = [8.430737686777691, 6.7429037100800215, 5.507118373636474]
renderView1.CameraFocalPoint = [-1.7464481617227953e-16, 0.9982166950000001, -6.565594592942839e-17]
renderView1.CameraViewUp = [-0.385105062916472, 0.8673653126663391, -0.3152324616839518]
renderView1.CameraParallelScale = 3.0005949060439203

# save animation
SaveAnimation(str(path_folder)+'/growth_factor_t.avi', renderView1, ImageResolution=[2000, 1352],
    OverrideColorPalette='WhiteBackground',
    FrameRate=10,
    FrameWindow=[0, len(res)])

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'growth_factor', 'Y'))

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(growth_factorLUT, clip1Display)

# Rescale transfer function
growth_factorLUT.RescaleTransferFunction(0.9737437445549506, 1.23539)

# Rescale transfer function
growth_factorPWF.RescaleTransferFunction(0.9737437445549506, 1.23539)

# Rescale transfer function
growth_factorLUT.RescaleTransferFunction(1.0, 2.0)

# Rescale transfer function
growth_factorPWF.RescaleTransferFunction(1.0, 2.0)

# layout/tab size in pixels
layout1.SetSize(2000, 1354)

# current camera placement for renderView1
renderView1.CameraPosition = [8.430737686777691, 6.7429037100800215, 5.507118373636474]
renderView1.CameraFocalPoint = [-1.7464481617227953e-16, 0.9982166950000001, -6.565594592942839e-17]
renderView1.CameraViewUp = [-0.385105062916472, 0.8673653126663391, -0.3152324616839518]
renderView1.CameraParallelScale = 3.0005949060439203

# save animation
SaveAnimation(str(path_folder)+'/growth_factor_r.avi', renderView1, ImageResolution=[2000, 1352],
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
renderView1.CameraPosition = [8.430737686777691, 6.7429037100800215, 5.507118373636474]
renderView1.CameraFocalPoint = [-1.7464481617227953e-16, 0.9982166950000001, -6.565594592942839e-17]
renderView1.CameraViewUp = [-0.385105062916472, 0.8673653126663391, -0.3152324616839518]
renderView1.CameraParallelScale = 3.0005949060439203

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
