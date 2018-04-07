# state file generated using paraview version 5.3.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1007, 411]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [50.0, 250.0, 0.0]
renderView1.StereoType = 0
renderView1.CameraPosition = [28.08200175660954, 107.88073302294208, 10000.0]
renderView1.CameraFocalPoint = [28.08200175660954, 107.88073302294208, 0.0]
renderView1.CameraViewUp = [1.0, 2.220446049250313e-16, 0.0]
renderView1.CameraParallelScale = 37.896837993316716
renderView1.Background = [0.09411764705882353, 0.10196078431372549, 0.12941176470588237]

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
data00 = LegacyVTKReader(FileNames=['/home/ema/simulazioni/disch_outcap/out/data.0000.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0001.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0002.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0003.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0004.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0005.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0006.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0007.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0008.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0009.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0010.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0011.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0012.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0013.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0014.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0015.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0016.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0017.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0018.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0019.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0020.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0021.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0022.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0023.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0024.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0025.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0026.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0027.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0028.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0029.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0030.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0031.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0032.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0033.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0034.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0035.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0036.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0037.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0038.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0039.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0040.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0041.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0042.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0043.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0044.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0045.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0046.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0047.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0048.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0049.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0050.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0051.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0052.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0053.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0054.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0055.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0056.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0057.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0058.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0059.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0060.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0061.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0062.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0063.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0064.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0065.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0066.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0067.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0068.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0069.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0070.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0071.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0072.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0073.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0074.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0075.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0076.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0077.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0078.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0079.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0080.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0081.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0082.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0083.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0084.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0085.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0086.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0087.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0088.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0089.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0090.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0091.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0092.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0093.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0094.vtk', '/home/ema/simulazioni/disch_outcap/out/data.0095.vtk'])

# create a new 'Glyph'
glyph1 = Glyph(Input=data00,
    GlyphType='Sphere')
glyph1.Scalars = ['CELLS', 'interBound']
glyph1.Vectors = ['CELLS', '2D_Velocity_Field']
glyph1.ScaleFactor = 2.0
glyph1.GlyphMode = 'All Points'
glyph1.GlyphTransform = 'Transform2'

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'bx3'
bx3LUT = GetColorTransferFunction('bx3')
bx3LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 100.0, 0.865003, 0.865003, 0.865003, 200.0, 0.705882, 0.0156863, 0.14902]
bx3LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'bx3'
bx3PWF = GetOpacityTransferFunction('bx3')
bx3PWF.Points = [0.0, 0.0, 0.5, 0.0, 200.0, 1.0, 0.5, 0.0]
bx3PWF.ScalarRangeInitialized = 1

# get color transfer function/color map for 'interBound'
interBoundLUT = GetColorTransferFunction('interBound')
interBoundLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 16.649999618530273, 0.865003, 0.865003, 0.865003, 33.29999923706055, 0.705882, 0.0156863, 0.14902]
interBoundLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'interBound'
interBoundPWF = GetOpacityTransferFunction('interBound')
interBoundPWF.Points = [0.0, 0.0, 0.5, 0.0, 33.29999923706055, 1.0, 0.5, 0.0]
interBoundPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from data00
data00Display = Show(data00, renderView1)
# trace defaults for the display properties.
data00Display.Representation = 'Surface With Edges'
data00Display.ColorArrayName = ['CELLS', 'bx3']
data00Display.LookupTable = bx3LUT
data00Display.EdgeColor = [0.7529411764705882, 0.7647058823529411, 0.6862745098039216]
data00Display.OSPRayScaleArray = 'rho'
data00Display.OSPRayScaleFunction = 'PiecewiseFunction'
data00Display.SelectOrientationVectors = '2D_Velocity_Field'
data00Display.ScaleFactor = 50.0
data00Display.SelectScaleArray = 'rho'
data00Display.GlyphType = 'Arrow'
data00Display.PolarAxes = 'PolarAxesRepresentation'
data00Display.GaussianRadius = 25.0
data00Display.SetScaleArray = [None, '']
data00Display.ScaleTransferFunction = 'PiecewiseFunction'
data00Display.OpacityArray = [None, '']
data00Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
data00Display.SetScalarBarVisibility(renderView1, True)

# show data from glyph1
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = ['POINTS', 'interBound']
glyph1Display.LookupTable = interBoundLUT
glyph1Display.EdgeColor = [0.7529411764705882, 0.7647058823529411, 0.6862745098039216]
glyph1Display.OSPRayScaleArray = 'GlyphScale'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'GlyphVector'
glyph1Display.ScaleFactor = 50.29577552080155
glyph1Display.SelectScaleArray = 'GlyphScale'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
glyph1Display.GaussianRadius = 25.147887760400774
glyph1Display.SetScaleArray = ['POINTS', 'GlyphScale']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'GlyphScale']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# setup the color legend parameters for each legend in this view

# get color legend/bar for interBoundLUT in view renderView1
interBoundLUTColorBar = GetScalarBar(interBoundLUT, renderView1)
interBoundLUTColorBar.Title = 'interBound'
interBoundLUTColorBar.ComponentTitle = ''

# get color legend/bar for bx3LUT in view renderView1
bx3LUTColorBar = GetScalarBar(bx3LUT, renderView1)
bx3LUTColorBar.Position = [0.85, 0.52]
bx3LUTColorBar.Position2 = [0.12, 0.42999999999999994]
bx3LUTColorBar.Title = 'bx3'
bx3LUTColorBar.ComponentTitle = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(data00)
# ----------------------------------------------------------------