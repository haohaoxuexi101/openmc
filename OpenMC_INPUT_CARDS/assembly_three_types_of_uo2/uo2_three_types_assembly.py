from math import pi

import openmc
import openmc.deplete
import numpy as np



uo2 = openmc.Material(name='UO2 fuel1')
uo2.temperature = 900.0
uo2.set_density('g/cm3', 10.412)

uo2.add_element('U', 1, enrichment=2, enrichment_type='wo')
uo2.add_element('O', 2)

xuo2 = openmc.Material(name='UO2 fuel2')
xuo2.temperature = 900.0
xuo2.set_density('g/cm3', 10.412)

xuo2.add_element('U', 1, enrichment=4, enrichment_type='wo')
xuo2.add_element('O', 2)

yuo2 = openmc.Material(name='UO2 fuel3')
yuo2.temperature = 900.0
yuo2.set_density('g/cm3', 10.412)

yuo2.add_element('U', 1, enrichment=6, enrichment_type='wo')
yuo2.add_element('O', 2)

helium = openmc.Material(name='Helium for gap')
helium.temperature = 800
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)

ss304 = openmc.Material(name='SS304')
ss304.temperature = 560
ss304.set_density('g/cc', 6.6)
ss304.add_element('Zr', 1.0)

borated_water = openmc.model.borated_water(680.916772554002,temperature=594.14,temp_unit='K',density=0.6775, name = "water")
borated_water.temperature=594.14
borated_water.add_s_alpha_beta('c_H_in_H2O')

# materials = openmc.Materials([gd_o_uo2, gd_uo2, uo2, helium, ss304, borated_water])
# materials.export_to_xml()
###############################################################################
#                             Create geometry
###############################################################################

# Define surfaces
pin_pitch =1.26
fuel_inner_radius = 0.0
fuel_outer_radius =0.805/2
clad_inner_radius = 0.822/2
clad_outer_radius = 0.95/2

fuel_or = openmc.ZCylinder(r=fuel_outer_radius, name='Fuel OR')
clad_ir = openmc.ZCylinder(r=clad_inner_radius, name='Clad IR')
clad_or = openmc.ZCylinder(r=clad_outer_radius, name='Clad OR')

gt_ir = openmc.ZCylinder(r=0.5715, name='guide tube IR')
gt_or = openmc.ZCylinder(r=0.6121, name='guide tube OR')


# Define uo2 cells
fuel = openmc.Cell(fill=uo2, region=-fuel_or, name = 'UO2 fuel1 cell')
gap = openmc.Cell(fill=helium, region=+fuel_or & -clad_ir, name = "gap1")
clad = openmc.Cell(fill=ss304, region=+clad_ir & -clad_or, name = "clad1")
water = openmc.Cell(fill=borated_water, region=+clad_or, name = "coolant1")

# Define uo2 cells
fuel1 = openmc.Cell(fill=xuo2, region=-fuel_or, name = 'UO2 fuel2 cell')
gap1 = openmc.Cell(fill=helium, region=+fuel_or & -clad_ir, name = "gap2")
clad1 = openmc.Cell(fill=ss304, region=+clad_ir & -clad_or, name = "clad2")
water1 = openmc.Cell(fill=borated_water, region=+clad_or, name = "coolant2")

# Define uo2 cells
fuel2 = openmc.Cell(fill=yuo2, region=-fuel_or, name = 'UO2 fuel3 cell')
gap2 = openmc.Cell(fill=helium, region=+fuel_or & -clad_ir, name = "gap3")
clad2 = openmc.Cell(fill=ss304, region=+clad_ir & -clad_or, name = "clad3")
water2 = openmc.Cell(fill=borated_water, region=+clad_or, name = "coolant3")


materials = openmc.Materials([uo2, xuo2, yuo2, helium, ss304, borated_water])
materials.export_to_xml()

# gd_o_gap = openmc.Cell(fill=helium, region=+fuel_or & -clad_ir)
# gd_o_clad = openmc.Cell(fill=ss304, region=+clad_ir & -clad_or)
# gd_o_water = openmc.Cell(fill=borated_water, region=+clad_or)


# Define pin universe
uo2_universe = openmc.Universe(cells=[fuel, gap, clad, water], name = "uo2_universe")
xuo2_universe = openmc.Universe(cells=[fuel1, gap1, clad1, water1], name = "xuo2_universe")
yuo2_universe = openmc.Universe(cells=[fuel2, gap2, clad2, water2], name = "yuo2_universe")


pin_boundary = openmc.model.RectangularPrism(pin_pitch, pin_pitch)
assembly_pitch = 21.42
# assembly_pitch = 21.4
# 这里如果不更改，会造成组件最右侧的燃料栅元的最右侧边界位置并不是pitch/2，而是pitch/2 - 0.02
# 问题是为什么最左侧没事？（负值的情况）
# lattice的设置规则是左下角的位置是lattice的lower_left属性，而不是center属性？
assembly_boundary = openmc.model.RectangularPrism(assembly_pitch, assembly_pitch, boundary_type='reflective')

fuel_lat = openmc.RectLattice(name="sf96_1 assembly")
fuel_lat.center = (0., 0.)
fuel_lat.pitch = (pin_pitch, pin_pitch)
fuel_lat.lower_left = [-assembly_pitch/2, -assembly_pitch/2]
# coolant universe for empty pin space filling
coolant_cell = openmc.Cell(fill=borated_water, name = "lattice outer water cell")
coolant_u = openmc.Universe(cells=[coolant_cell], name = "lattice outer water universe")
fuel_lat.outer = coolant_u

u0 = uo2_universe
ux = xuo2_universe
uy = yuo2_universe

fuel_lat.universes = [  
                        [u0, u0, u0, u0, u0, u0, u0, u0, u0, uy, uy, uy, uy, uy, uy, uy, uy],
                        [u0, u0, u0, u0, u0, u0, u0, u0, u0, uy, uy, uy, uy, uy, uy, uy, uy],
                        [u0, u0, u0, u0, u0, u0, u0, u0, u0, uy, uy, uy, uy, uy, uy, uy, uy],
                        [u0, u0, u0, u0, u0, u0, u0, u0, u0, uy, uy, uy, uy, uy, uy, uy, uy],
                        [u0, u0, u0, u0, u0, u0, u0, u0, u0, uy, uy, uy, uy, uy, uy, uy, uy],
                        [u0, u0, u0, u0, u0, u0, u0, u0, u0, uy, uy, uy, uy, uy, uy, uy, uy],
                        [u0, u0, u0, u0, u0, u0, u0, u0, u0, uy, uy, uy, uy, uy, uy, uy, uy],
                        [u0, u0, u0, u0, u0, u0, u0, u0, u0, uy, uy, uy, uy, uy, uy, uy, uy],
                        [u0, u0, u0, u0, u0, u0, u0, u0, u0, u0, u0, u0, u0, u0, u0, u0, u0],
                        [ux, ux, ux, ux, ux, ux, ux, ux, u0, u0, u0, u0, u0, u0, u0, u0, u0],
                        [ux, ux, ux, ux, ux, ux, ux, ux, u0, u0, u0, u0, u0, u0, u0, u0, u0],
                        [ux, ux, ux, ux, ux, ux, ux, ux, u0, u0, u0, u0, u0, u0, u0, u0, u0],
                        [ux, ux, ux, ux, ux, ux, ux, ux, u0, u0, u0, u0, u0, u0, u0, u0, u0],
                        [ux, ux, ux, ux, ux, ux, ux, ux, u0, u0, u0, u0, u0, u0, u0, u0, u0],
                        [ux, ux, ux, ux, ux, ux, ux, ux, u0, u0, u0, u0, u0, u0, u0, u0, u0],
                        [ux, ux, ux, ux, ux, ux, ux, ux, u0, u0, u0, u0, u0, u0, u0, u0, u0],
                        [ux, ux, ux, ux, ux, ux, ux, ux, u0, u0, u0, u0, u0, u0, u0, u0, u0],
                    ]

fuel_lat_cell = openmc.Cell(fill=fuel_lat,region=-assembly_boundary, name = "fuel lattice")

geometry = openmc.Geometry([fuel_lat_cell])
geometry.export_to_xml()
fuel_lat_cell.plot(origin=(0., 0., 0.), pixels=(1500, 1500), color_by='cell')


SHEM_361 = np.array([
 1.1000E-10, 2.4999E-09, 4.5560E-09, 7.1453E-09, 1.0451E-08, 1.4830E-08, 2.0010E-08, 2.4939E-08, 2.9299E-08, 3.4400E-08, 4.0300E-08, 4.7302E-08,
 5.5498E-08, 6.5199E-08, 7.6497E-08, 8.9797E-08, 1.0430E-07, 1.1999E-07, 1.3800E-07, 1.6190E-07, 1.9000E-07, 2.0961E-07, 2.3119E-07, 2.5500E-07,
 2.7999E-07, 3.0501E-07, 3.2501E-07, 3.5299E-07, 3.9000E-07, 4.3158E-07, 4.7502E-07, 5.2001E-07, 5.5499E-07, 5.9499E-07, 6.2500E-07, 7.2000E-07,
 8.2004E-07, 8.8002E-07, 9.1998E-07, 9.4402E-07, 9.6396E-07, 9.8196E-07, 9.9650E-07, 1.0090E-06, 1.0210E-06, 1.0350E-06, 1.0780E-06, 1.0920E-06,
 1.1040E-06, 1.1160E-06, 1.1300E-06, 1.1480E-06, 1.1700E-06, 1.2140E-06, 1.2509E-06, 1.2930E-06, 1.3310E-06, 1.3810E-06, 1.4100E-06, 1.4440E-06,
 1.5200E-06, 1.5880E-06, 1.6689E-06, 1.7800E-06, 1.9001E-06, 1.9899E-06, 2.0701E-06, 2.1569E-06, 2.2171E-06, 2.2730E-06, 2.3301E-06, 2.4699E-06,
 2.5500E-06, 2.5901E-06, 2.6201E-06, 2.6400E-06, 2.7001E-06, 2.7199E-06, 2.7409E-06, 2.7751E-06, 2.8840E-06, 3.1421E-06, 3.5431E-06, 3.7121E-06,
 3.8822E-06, 4.0000E-06, 4.2198E-06, 4.3098E-06, 4.4198E-06, 4.7678E-06, 4.9332E-06, 5.1100E-06, 5.2101E-06, 5.3201E-06, 5.3800E-06, 5.4102E-06,
 5.4882E-06, 5.5300E-06, 5.6198E-06, 5.7201E-06, 5.8002E-06, 5.9601E-06, 6.0599E-06, 6.1601E-06, 6.2802E-06, 6.3598E-06, 6.4321E-06, 6.4818E-06,
 6.5149E-06, 6.5391E-06, 6.5561E-06, 6.5718E-06, 6.5883E-06, 6.6061E-06, 6.6313E-06, 6.7167E-06, 6.7423E-06, 6.7598E-06, 6.7761E-06, 6.7917E-06,
 6.8107E-06, 6.8353E-06, 6.8702E-06, 6.9178E-06, 6.9943E-06, 7.1399E-06, 7.3802E-06, 7.6004E-06, 7.7399E-06, 7.8397E-06, 7.9701E-06, 8.1303E-06,
 8.3003E-06, 8.5241E-06, 8.6737E-06, 8.8004E-06, 8.9800E-06, 9.1403E-06, 9.5000E-06, 1.0579E-05, 1.0804E-05, 1.1053E-05, 1.1269E-05, 1.1589E-05,
 1.1709E-05, 1.1815E-05, 1.1979E-05, 1.2130E-05, 1.2309E-05, 1.2472E-05, 1.2600E-05, 1.3330E-05, 1.3546E-05, 1.4050E-05, 1.4251E-05, 1.4470E-05,
 1.4595E-05, 1.4730E-05, 1.4866E-05, 1.5779E-05, 1.6050E-05, 1.6550E-05, 1.6831E-05, 1.7446E-05, 1.7565E-05, 1.7759E-05, 1.7959E-05, 1.9085E-05,
 1.9200E-05, 1.9393E-05, 1.9597E-05, 2.0073E-05, 2.0275E-05, 2.0418E-05, 2.0520E-05, 2.0602E-05, 2.0685E-05, 2.0768E-05, 2.0976E-05, 2.1060E-05,
 2.1145E-05, 2.1230E-05, 2.1336E-05, 2.1486E-05, 2.1702E-05, 2.2001E-05, 2.2156E-05, 2.2378E-05, 2.2536E-05, 2.4658E-05, 2.7885E-05, 3.1693E-05,
 3.3085E-05, 3.4539E-05, 3.5698E-05, 3.6057E-05, 3.6419E-05, 3.6859E-05, 3.7304E-05, 3.7792E-05, 3.8787E-05, 3.9730E-05, 4.1227E-05, 4.2144E-05,
 4.3125E-05, 4.4172E-05, 4.5290E-05, 4.6205E-05, 4.7517E-05, 4.9259E-05, 5.1785E-05, 5.2990E-05, 5.4060E-05, 5.7059E-05, 5.9925E-05, 6.2308E-05,
 6.3631E-05, 6.4592E-05, 6.5046E-05, 6.5503E-05, 6.5831E-05, 6.6161E-05, 6.6493E-05, 6.6826E-05, 6.9068E-05, 7.1887E-05, 7.3559E-05, 7.6332E-05,
 7.9368E-05, 8.3939E-05, 8.8774E-05, 9.3326E-05, 9.7329E-05, 1.0059E-04, 1.0110E-04, 1.0161E-04, 1.0211E-04, 1.0304E-04, 1.0565E-04, 1.1029E-04,
 1.1285E-04, 1.1548E-04, 1.1652E-04, 1.1758E-04, 1.2055E-04, 1.2623E-04, 1.3270E-04, 1.3950E-04, 1.4666E-04, 1.5418E-04, 1.6306E-04, 1.6752E-04,
 1.7523E-04, 1.8329E-04, 1.8495E-04, 1.8625E-04, 1.8756E-04, 1.8888E-04, 1.9020E-04, 1.9308E-04, 1.9600E-04, 2.0096E-04, 2.1211E-04, 2.2432E-04,
 2.3559E-04, 2.4180E-04, 2.5675E-04, 2.6830E-04, 2.7647E-04, 2.8489E-04, 2.8833E-04, 2.9592E-04, 3.1993E-04, 3.3532E-04, 3.5357E-04, 3.7170E-04,
 3.9076E-04, 4.1909E-04, 4.5400E-04, 5.0175E-04, 5.3920E-04, 5.7715E-04, 5.9294E-04, 6.0010E-04, 6.1283E-04, 6.4684E-04, 6.7729E-04, 7.4852E-04,
 8.3222E-04, 9.0968E-04, 9.8249E-04, 1.0643E-03, 1.1347E-03, 1.3436E-03, 1.5862E-03, 1.8118E-03, 2.0841E-03, 2.3973E-03, 2.7002E-03, 2.9962E-03,
 3.4811E-03, 4.0973E-03, 5.0045E-03, 6.1125E-03, 7.4658E-03, 9.1188E-03, 1.1138E-02, 1.3604E-02, 1.4900E-02, 1.6200E-02, 1.8585E-02, 2.2699E-02,
 2.4999E-02, 2.6100E-02, 2.7394E-02, 2.9281E-02, 3.3460E-02, 3.6979E-02, 4.0868E-02, 4.9916E-02, 5.5166E-02, 6.7379E-02, 8.2297E-02, 9.4665E-02,
 1.1562E-01, 1.2277E-01, 1.4010E-01, 1.6507E-01, 1.9507E-01, 2.3006E-01, 2.6783E-01, 3.2065E-01, 3.8388E-01, 4.1250E-01, 4.5602E-01, 4.9400E-01,
 5.7844E-01, 7.0651E-01, 8.6001E-01, 9.5112E-01, 1.0511E+00, 1.1620E+00, 1.2870E+00, 1.3369E+00, 1.4058E+00, 1.6365E+00, 1.9014E+00, 2.2313E+00,
 2.7253E+00, 3.3287E+00, 4.0657E+00, 4.9658E+00, 6.0653E+00, 6.7032E+00, 7.4082E+00, 8.1873E+00, 9.0484E+00, 1.0000E+01, 1.1618E+01, 1.3840E+01,
 1.4918E+01, 1.9640E+01])*1E6


whole_geometry_filter = openmc.CellFilter(fuel_lat_cell)
fuel_filter = openmc.MaterialFilter([uo2, xuo2, yuo2])

## definition for xslib tallies

# whole power tally for xs
overall_power_tally = openmc.Tally(name="_tally_overall_power", tally_id=40000)
overall_power_tally.scores = ["heating-local", "fission-q-recoverable", "fission-q-prompt", "kappa-fission", "fission", '(n,gamma)']
overall_power_tally.filters = [whole_geometry_filter]

fuel_tally = openmc.Tally()
fuel_tally.scores = ["heating-local", "fission-q-recoverable", "fission-q-prompt", "kappa-fission", "fission", '(n,gamma)']
fuel_tally.filters = [fuel_filter]

tallies = openmc.Tallies([overall_power_tally, fuel_tally])
tallies.export_to_xml()

###############################################################################
#                     Set volumes of depletable materials
###############################################################################

# Set material volume for depletion. For 2D simulations, this should be an area.
uo2.volume = (289 - 2 * 8 * 8) * pi * fuel_or.r**2
xuo2.volume = 8 * 8 * pi * fuel_or.r**2
yuo2.volume = 8 * 8 * pi * fuel_or.r**2

###############################################################################
#                     Transport calculation settings
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings = openmc.Settings()
settings.batches = 20
settings.inactive = 5
settings.particles = 10000

# settings.photon_transport = True

settings.temperature = {'tolerance': 1000}

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-assembly_pitch/2, -assembly_pitch/2, -1, assembly_pitch/2, assembly_pitch/2, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.source.IndependentSource(space=uniform_dist, constraints={'fissionable': True})
settings.export_to_xml()
###############################################################################
#                   Initialize and run depletion calculation
###############################################################################

# Create depletion "operator"
chain_file = './chain_casl_pwr.xml'
model = openmc.Model(geometry=geometry, settings=settings)
# op = openmc.deplete.CoupledOperator(model, chain_file,
#                              normalization_mode="fission-q", fission_q=serpent_fq,
#                              fission_yield_mode="average")
op = openmc.deplete.CoupledOperator(model, chain_file,
                             normalization_mode="energy-deposition",
                             fission_yield_mode="average")

# Perform simulation using the predictor algorithm
time_steps = [20]*5
power_density = [20]*5
integrator = openmc.deplete.PredictorIntegrator(op, time_steps, power_density=power_density, timestep_units='d')
# integrator = openmc.deplete.CECMIntegrator(op, time_steps, power_density=power_density, timestep_units='d')
integrator.integrate()
# openmc.run()
