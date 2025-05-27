import openmc
from matplotlib.colors import LogNorm
from math import pi

from openmc_regular_mesh_plotter import plot_mesh_tally

# MATERIALS
mat_1 = openmc.Material()
mat_1.add_element("N", 0.7)
mat_1.add_element("O", 0.3)
mat_1.set_density("kg/m3", 1.29)

mat_2 = openmc.Material()
mat_2.add_element("N", 0.7)
mat_2.add_element("O", 0.3)
mat_2.set_density("kg/m3", 1.29)
my_materials = openmc.Materials([mat_1, mat_2])

# GEOMETRY
# surfaces
inner_surface = openmc.Sphere(r=300)
outer_surface = openmc.Sphere(r=400, boundary_type="vacuum", surface_id = 3)
# regions
inner_region = -inner_surface
outer_region = -outer_surface & +inner_surface
# cells
inner_cell = openmc.Cell(region=inner_region)
inner_cell.fill = mat_2
outer_cell = openmc.Cell(region=outer_region)
outer_cell.fill = mat_1
my_geometry = openmc.Geometry([inner_cell, outer_cell])

# SIMULATION SETTINGS
my_settings = openmc.Settings()
my_settings.batches = 20
my_settings.inactive = 0
my_settings.particles = 500000
# my_settings.photon_transport = True
my_settings.run_mode = "fixed source"
my_settings.create_fission_neutrons = False
# Create a DT point source
try:
    source = openmc.IndependentSource()
except:
    # work with older versions of openmc
    source = openmc.Source()
source.space = openmc.stats.Point((-400, 0, 0))
# source.angle = openmc.stats.Isotropic()
phi = openmc.stats.Uniform(0.0*pi, 2.0*pi)
mu = openmc.stats.Uniform(0.7, 1.0)
angle = openmc.stats.PolarAzimuthal(mu=mu, phi=phi, reference_uvw=(1.0, 0.0, 0.0))
source.angle = angle
source.energy = openmc.stats.Discrete([1e4], [1])
source.particle = "neutron"
my_settings.source = source

my_settings.surf_source_write = {
    'surface_ids': [3],
    'max_particles': 5000000
    }

# Tallies
my_tallies = openmc.Tallies()
mesh = openmc.RegularMesh().from_domain(
    my_geometry,  # the corners of the mesh are being set automatically to surround the geometry
    dimension=[40, 40, 1],
)
mesh_filter = openmc.MeshFilter(mesh)
mesh_tally_1 = openmc.Tally(name="mesh_tally")
mesh_tally_1.filters = [mesh_filter]
mesh_tally_1.scores = ["heating"]
my_tallies.append(mesh_tally_1)

model = openmc.model.Model(my_geometry, my_materials, my_settings, my_tallies)
sp_filename = model.run()

# post process simulation result
statepoint = openmc.StatePoint(sp_filename)

# extracts the mesh tally by name
my_mesh_tally = statepoint.get_tally(name="mesh_tally")

# default tally units for heating are in eV per source neutron
# for this example plot we want Mega Joules per second per cm3 or Mjcm^-3s^-1
neutrons_per_second = 1e21
eV_to_joules = 1.60218e-19
joules_to_mega_joules = 1e-6
scaling_factor = neutrons_per_second * eV_to_joules * joules_to_mega_joules
# note that volume_normalization is enabled so this will also change the units to divide by the volume of each mesh voxel
# alternatively you could set volume_normalization to false and divide by the mesh.volume[0][0][0] in the scaling factor
# in a regular mesh all the voxels have the same volume so the [0][0][0] just picks the first volume

plot = plot_mesh_tally(
    basis="xy",  # as the mesh dimention is [1,40,40] only the yz basis can be plotted
    tally=my_mesh_tally,
    # outline=True,  # enables an outline around the geometry
    geometry=my_geometry,  # needed for outline
    norm=LogNorm(),  # log scale
    colorbar=False,
)

plot.figure.savefig("example_openmc_2d_regular_mesh_plotter.png")
plot.title.set_text("")