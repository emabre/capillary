# state file generated using paraview version 5.4.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1372, 723]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [30.0, 60.0, 0.0]
renderView1.StereoType = 0
renderView1.CameraPosition = [32.949863117104435, 97.50263715651157, -259.185097056601]
renderView1.CameraFocalPoint = [32.949863117104435, 97.50263715651157, 0.0]
renderView1.CameraViewUp = [4.44089209850063e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 118.83992466862519
renderView1.Background = [0.1411764705882353, 0.1411764705882353, 0.17647058823529413]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
datavtk = LegacyVTKReader(FileNames=['/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0000.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0001.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0002.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0003.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0004.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0005.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0006.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0007.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0008.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0009.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0010.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0011.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0012.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0013.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0014.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0015.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0016.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0017.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0018.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0019.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0020.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0021.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0022.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0023.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0024.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0025.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0026.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0027.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0028.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0029.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0030.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0031.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0032.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0033.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0034.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0035.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0036.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0037.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0038.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0039.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0040.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0041.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0042.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0043.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0044.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0045.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0046.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0047.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0048.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0049.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0050.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0051.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0052.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0053.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0054.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0055.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0056.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0057.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0058.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0059.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0060.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0061.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0062.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0063.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0064.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0065.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0066.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0067.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0068.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0069.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0070.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0071.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0072.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0073.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0074.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0075.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0076.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0077.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0078.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0079.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0080.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0081.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0082.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0083.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0084.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0085.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0086.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0087.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0088.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0089.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0090.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0091.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0092.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0093.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0094.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0095.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0096.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0097.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0098.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0099.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0100.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0101.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0102.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0103.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0104.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0105.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0106.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0107.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0108.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0109.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0110.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0111.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0112.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0113.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0114.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0115.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0116.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0117.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0118.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0119.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0120.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0121.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0122.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0123.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0124.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0125.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0126.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0127.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0128.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0129.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0130.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0131.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0132.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0133.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0134.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0135.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0136.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0137.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0138.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0139.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0140.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0141.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0142.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0143.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0144.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0145.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0146.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0147.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0148.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0149.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0150.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0151.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0152.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0153.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0154.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0155.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0156.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0157.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0158.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0159.vtk', '/home/konrad/simulazioni/sims_pluto/disch_outcap/out/data.0160.vtk'])

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'T'
tLUT = GetColorTransferFunction('T')
tLUT.RGBPoints = [1202.72204589844, 0.231373, 0.298039, 0.752941, 2901.36102294922, 0.865003, 0.865003, 0.865003, 4600.0, 0.705882, 0.0156863, 0.14902]
tLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'T'
tPWF = GetOpacityTransferFunction('T')
tPWF.Points = [1202.72204589844, 0.0, 0.5, 0.0, 4600.0, 1.0, 0.5, 0.0]
tPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from datavtk
datavtkDisplay = Show(datavtk, renderView1)
# trace defaults for the display properties.
datavtkDisplay.Representation = 'Surface With Edges'
datavtkDisplay.ColorArrayName = ['CELLS', 'T']
datavtkDisplay.LookupTable = tLUT
datavtkDisplay.EdgeColor = [0.7686274509803922, 0.7843137254901961, 0.6941176470588235]
datavtkDisplay.Scale = [1.0, 0.1, 1.0]
datavtkDisplay.OSPRayScaleArray = 'rho'
datavtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
datavtkDisplay.SelectOrientationVectors = 'None'
datavtkDisplay.ScaleFactor = 120.0
datavtkDisplay.SelectScaleArray = 'rho'
datavtkDisplay.GlyphType = 'Arrow'
datavtkDisplay.GlyphTableIndexArray = 'rho'
datavtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
datavtkDisplay.PolarAxes = 'PolarAxesRepresentation'
datavtkDisplay.GaussianRadius = 60.0
datavtkDisplay.SetScaleArray = [None, '']
datavtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
datavtkDisplay.OpacityArray = [None, '']
datavtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
datavtkDisplay.OSPRayScaleFunction.Points = [5.32782712934216e-15, 0.0, 0.5, 0.0, 1e+30, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
datavtkDisplay.PolarAxes.Scale = [1.0, 0.1, 1.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
datavtkDisplay.ScaleTransferFunction.Points = [5.32782712934216e-15, 0.0, 0.5, 0.0, 1e+30, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
datavtkDisplay.OpacityTransferFunction.Points = [5.32782712934216e-15, 0.0, 0.5, 0.0, 1e+30, 1.0, 0.5, 0.0]

# show color legend
datavtkDisplay.SetScalarBarVisibility(renderView1, True)

# setup the color legend parameters for each legend in this view

# get color legend/bar for tLUT in view renderView1
tLUTColorBar = GetScalarBar(tLUT, renderView1)
tLUTColorBar.Title = 'T'
tLUTColorBar.ComponentTitle = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------