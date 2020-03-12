#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import numpy as np
from six.moves import cPickle as pickle
from math import pi

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
# homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
udir = 'TIDY/ARCHIVE/201[34567]/d05/U/'
vdir = 'TIDY/ARCHIVE/201[34567]/d05/V/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# --------------------------------------------------------------------------- #

# Load the grid spacings.
e2u = nemo.load_field('e2t', homedir, nemodir, gridfile, 'U')

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T')[:, :, [32, 52]])
umask3d = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))
umask2d = umask3d[:, :, [32, 52]]
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V')[:, :, [32, 52]])
fmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'Z')[:, :, 0])

# --------------------------------------------------------------------------- #
# Find the files that I want to know the turning angle for.

ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*'])))
vfiles = sorted(glob.glob(''.join([homedir, nemodir, vdir, '*'])))

# --------------------------------------------------------------------------- #

# Preallocate the output dictionary.
delta_angle = { 'boundaries' : np.zeros(shape=[len(ufiles), 2]),
                'domain_037' : np.zeros(shape=[37, len(ufiles)]),
                'north_037'  : np.zeros(shape=[37, len(ufiles)]),
                'acc_037'    : np.zeros(shape=[37, len(ufiles)]),
                'south_037'  : np.zeros(shape=[37, len(ufiles)]),
                'domain_361' : np.zeros(shape=[361, len(ufiles)]),
                'north_361'  : np.zeros(shape=[361, len(ufiles)]),
                'acc_361'    : np.zeros(shape=[361, len(ufiles)]),
                'south_361'  : np.zeros(shape=[361, len(ufiles)]) }

# Set the edges for the histogram.
edges = np.linspace(-185., 185., 38)

# Loop over the U/V files and calculate the turning angle..
for k in range(len(ufiles)):
    print(ufiles[k])
 
    # Load the U velocity field and remask it.
    uVel = np.squeeze(nemo.load_field('vozoce3u', '', '', ufiles[k], 'U'))
    uVel = nemo.mask_field(uVel, umask3d)

    # Calculate the transport streamfunction to create a map of where to calculate the turning angle.
    psi = nemo.calc_psi(uVel.data, e2u, np.ones((4321,2000,75)), fmask, umask3d, nx, ny)

    # Select the high and low boundaries as specfied by psi at Drake Passage.
    delta_angle['boundaries'][k, :] = [ psi[2599,:].max() - 20.E6, psi[2599,:].min() + 20.E6 ]

    # North of the ACC.
    north = psi > delta_angle['boundaries'][k,0]

    # ACC region.
    acc = (psi > delta_angle['boundaries'][k,1]) & (psi < delta_angle['boundaries'][k,0])

    # South of the ACC.
    south = psi < delta_angle['boundaries'][k,1]

    # Load the e3u for the two specific levels.
    e3u = np.squeeze(nemo.load_field('e3u', '', '', ufiles[k], 'U')[:,:,[32,52]])

    # Discard the unwanted bits of uVel. Factor out e3u & discard.
    uVel = uVel[:,:,[32,52]]
    uVel = uVel / e3u
    e3u = None

    # Load vVel, factor out e3v and discard.
    vVel = np.squeeze(nemo.load_field('vomece3v', '', '', vfiles[k], 'V')[:,:,[32,52]])
    e3v = np.squeeze(nemo.load_field('e3v', '', '', vfiles[k], 'V')[:,:,[32,52]])
    vVel = nemo.mask_field(vVel / e3v, vmask)
    e3v = None

    # Calculate uVel & vVel on T points.
    uVel = nemo.calc_u_on_t(uVel.data * umask2d, nx, ny, 2, tmask)
    vVel = nemo.calc_v_on_t(vVel.data * vmask, nx, ny, 2, tmask)

    # Take the dot_product of the velocities vectors at the two different levels.
    dot_product = uVel.data[:,:,0]*uVel.data[:,:,1] + vVel.data[:,:,0]*vVel.data[:,:,1]

    # Magnitude of the two vectors.
    mag_a = np.sqrt(uVel.data[:,:,0]*uVel.data[:,:,0] + vVel.data[:,:,0]*vVel.data[:,:,0])
    mag_b = np.sqrt(uVel.data[:,:,1]*uVel.data[:,:,1] + vVel.data[:,:,1]*vVel.data[:,:,1])

    # Cosine of the angle between the two vectors - cap magnitude at +/-1, then set to absolute value to prevent angles > 180.
    cos_theta = dot_product / ( mag_a * mag_b )
    cos_theta[ cos_theta > 1.0 ] = 1.
    cos_theta[ cos_theta < -1.0 ] = -1.
    cos_theta = np.absolute(cos_theta)

    # Angle between the two vectors.
    mag_angle = np.arccos( cos_theta )

    # Direction of angle (clockwise/anticlockwise).
    direction = np.sign( np.squeeze( uVel.data[:,:,0]*vVel.data[:,:,1] - vVel.data[:,:,0]*uVel.data[:,:,1] ) );
    
    # Save the angle, including the direction of rotation, for posterity.
    angle = mag_angle * direction
    
    # Set the edges for the histogram.
    edges = np.linspace(-185., 185., 38)

    # Save histograms in 10o bins.
    delta_angle['domain_037'][:,k] = np.histogram(180.0*angle[:]/pi, bins=edges)[0]
    delta_angle['north_037'][:,k] = np.histogram(180.0*angle[north[0:nx,0:ny]]/pi, bins=edges)[0]
    delta_angle['acc_037'][:,k] = np.histogram(180.0*angle[acc[0:nx,0:ny]]/pi, bins=edges)[0]
    delta_angle['south_037'][:,k] = np.histogram(180.0*angle[south[0:nx,0:ny]]/pi, bins=edges)[0]
    
    # Set the edges for the histogram.
    edges = np.linspace(-181., 181., 362)

    # Save histograms in 10o bins.
    delta_angle['domain_361'][:,k] = np.histogram(180.0*angle[:]/pi, bins=edges)[0]
    delta_angle['north_361'][:,k] = np.histogram(180.0*angle[north[0:nx,0:ny]]/pi, bins=edges)[0]
    delta_angle['acc_361'][:,k] = np.histogram(180.0*angle[acc[0:nx,0:ny]]/pi, bins=edges)[0]
    delta_angle['south_361'][:,k] = np.histogram(180.0*angle[south[0:nx,0:ny]]/pi, bins=edges)[0]

# --------------------------------------------------------------------------- #
# Spit the dictionary out to file, if requested.
# np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/delta_angle']), delta_angle)

with open(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/delta_angle.pkl']), 'wb') as f:
    pickle.dump(delta_angle, f)

# --------------------------------------------------------------------------- #
