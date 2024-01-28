from paraview.simple import *

# Load a PVD file
pvd_file = '/Users/tomoaki/SPH/Kamra2019_square_PS0d0075_CSML2d8_RK1_good/water_particle.pvd'
data_source = OpenDataFile(pvd_file)

# Create a render view
render_view = CreateView('RenderView')
render_view.ViewSize = [1920, 1080]  # You can adjust the size as needed

# Show data in the render view with PointGaussian representation
display = Show(data_source, render_view)
display.SetRepresentationType('PointGaussian')

# Adjust Point Gaussian properties (optional)
# For example, to change the size of the Gaussian blobs:
# display.GaussianRadius = 0.1

# Set up camera and other visualization parameters
# ...

# Create an animation
animation_scene = GetAnimationScene()
animation_scene.EndTime = 10  # Set the end time of your animation
animation_scene.Play()

# Save the animation
SaveAnimation('./animation.avi', render_view, FrameRate=15, FrameWindow=[0, 100])
