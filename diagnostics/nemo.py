#!/opt/local/bin/python
#
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# Define a function to average a T field to a V point.

def ave_t_onto_v(tfield, nx, ny, nz, vmask):

    # Import required commands.
    from numpy import ndarray

    # Define the output field as a numpy ndarray
    tfield_on_v = ndarray(shape=(nx, ny+1, nz))

    # Calculate the scalar product using 2nd order differences.
    tfield_on_v[0:nx, 1:ny, 0:nz] = 0.5*(tfield[0:nx, 1:ny, 0:nz] +
                                         tfield[0:nx, 0:ny-1, 0:nz])

    # Apply 'boundary' conditions at the north and south.
    tfield_on_v[0:nx, 0:1, 0:nz] = 0.0
    tfield_on_v[0:nx, -1:, 0:nz] = 0.0

    # Mask the output field.
    tfield_on_v = tfield_on_v * vmask

    return tfield_on_v


# --------------------------------------------------------------------------- #
# Define a function to average a U field to a V point.

def ave_u_onto_v(ufield, nx, ny, nz, vmask):

    # Import required commands.
    from numpy import ndarray

    # Define the output field as a numpy ndarray.
    ufield_on_v = ndarray(shape=(nx, ny+1, nz))

    # Calculate the scalar product using 2nd order differences.
    ufield_on_v[0:nx, 1:ny, 0:nz] = 0.25*(ufield[0:nx, 1:ny, 0:nz] +
                                          ufield[1:nx+1, 0:ny-1, 0:nz] +
                                          ufield[0:nx, 1:ny, 0:nz] +
                                          ufield[1:nx+1, 0:ny-1, 0:nz] )

    # Apply 'boundary' conditions at the north and south.
    ufield_on_v[0:nx, 0:1, 0:nz] = 0.0
    ufield_on_v[0:nx, -1:, 0:nz] = 0.0

    # Mask the output field.
    ufield_on_v = ufield_on_v * vmask

    return ufield_on_v


# --------------------------------------------------------------------------- #
# Define a function to average a U field to a V point.

def ave_v_onto_u(vfield, nx, ny, nz, umask):

    # Import required commands.
    from numpy import ndarray

    # Define the output field as a numpy ndarray.
    vfield_on_u = ndarray(shape=(nx+1, ny, nz))

    # Calculate the scalar product using 2nd order differences.
    vfield_on_u[1:nx, 0:ny, 0:nz] = 0.25*(vfield[1:nx, 0:ny, 0:nz] +
                                          vfield[1:nx, 1:ny+1, 0:nz] +
                                          vfield[0:nx-1, 0:ny, 0:nz] +
                                          vfield[0:nx-1, 1:ny+1, 0:nz] )

    # Apply 'boundary' conditions at the east and west.
    vfield_on_u[0, 0:ny, 0:nz] = 0.25*(vfield[0, 0:ny, 0:nz] +
                                       vfield[0, 1:ny+1, 0:nz] +
                                       vfield[nx-1, 0:ny, 0:nz] +
                                       vfield[nx-1, 1:ny+1, 0:nz] )
    vfield_on_u[nx, 0:ny, 0:nz] = vfield_on_u[0, 0:ny, 0:nz]

    # Mask the output field.
    vfield_on_u = vfield_on_u * umask

    return vfield_on_u


# --------------------------------------------------------------------------- #
# Define a function to construct an annual average from a specified
# set of years.

def calc_annual_ave(fieldname, fieldmask, fieldtype,
                    homedirectory, nemodirectory, filedir, prefix, suffix,
                    start_year=1948, no_years=1):

    # Import required modules.
    import glob

    # Initialise the output field to zero.
    annual = 0.0
    counter = 0.0

    # Loop over the number of years and print some shit to the screen
    for k in range(no_years):
        year = start_year + k

        # Find the names of the files.
        possible = ''.join([str(year), '/', filedir, prefix, str(year), '????', suffix])
        files = sorted(glob.glob(''.join([homedirectory, nemodirectory, possible])))

        # Load the current month & year combination into seasonal.
        for m in range(len(files)):
            print(files[m])
            annual = annual + load_field(fieldname, '',  '',
                                         files[m], fieldtype)
            counter = counter + 1

    # Divide by the number of files loaded.
    annual = annual / counter

    return annual


# --------------------------------------------------------------------------- #
# Define a function to construct an annual average from a specified
# set of years for the product of two terms.

def calc_annual_product(field1, field2, fieldmask, fieldtype,
                        homedirectory, nemodirectory, filedir, prefix, suffix,
                        start_year=1948, no_years=1):

    # Import required modules.
    import glob

    # Initialise the output field to zero.
    annual = 0.0
    counter = 0.0

    # Loop over the number of years and print some shit to the screen
    for k in range(no_years):
        year = start_year + k

        # Find the names of the files.
        possible = ''.join([str(year), '/', filedir, prefix, str(year),
                            '????', suffix])
        files = sorted(glob.glob(''.join([homedirectory, nemodirectory,
                                          possible])))

        # Load the current month & year combination into seasonal.
        for m in range(len(files)):
            print(files[m])
            annual = annual + \
                            (load_field(field1, '', '', files[m], fieldtype) *
                            load_field(field2, '', '', files[m], fieldtype))
            counter = counter + 1

    # Divide by the number of files loaded.
    annual = annual / counter

    return annual


# --------------------------------------------------------------------------- #
# Define a function to construct an annual average from a specified
# set of years for the ratio of two terms.

def calc_annual_ratio(field1, field2, fieldmask, fieldtype, exponent,
                      homedirectory, nemodirectory, filedir, prefix, suffix,
                      start_year=1948, no_years=1):

    # Import required modules.
    import glob
    from numpy import power

    # Initialise the output field to zero.
    annual = 0.0
    counter = 0.0

    # Loop over the number of years and print some shit to the screen
    for k in range(no_years):
        year = start_year + k

        # Find the names of the files.
        possible = ''.join([str(year), '/', filedir, prefix, str(year),
                            '????', suffix])
        files = sorted(glob.glob(''.join([homedirectory, nemodirectory,
                                          possible])))

        # Load the current month & year combination into seasonal.
        for m in range(len(files)):
            print(files[m])
            annual = annual + power \
                            (load_field(field1, '', '', files[m], fieldtype) /
                            load_field(field2, '', '', files[m], fieldtype),
                            exponent)
            counter = counter + 1

    # Divide by the number of files loaded.
    annual = annual / counter

    return annual


# --------------------------------------------------------------------------- #
# Define a function to calculate eddy kinetic energy from input velocities.

def calc_eke(uvelocity, uvel_sq, vvelocity, vvel_sq, nx, ny, nz, tmask):

    # Reynolds average the velocities.
    uu = uvel_sq - uvelocity * uvelocity
    vv = vvel_sq - vvelocity * vvelocity

    # Calculate eke using 2nd order differences.
    eke = 0.50 * calc_uv_on_t(uu, vv, nx, ny, nz, tmask)

    return eke


# --------------------------------------------------------------------------- #
# Define a function to calculate kinetic energy from input velocities.

def calc_ke(uvelocity, vvelocity, nx, ny, nz, tmask):

    # Square the velocities.
    uu = uvelocity * uvelocity
    vv = vvelocity * vvelocity

    # Calculate ke using 2nd order differences.
    ke = 0.50 * calc_uv_on_t(uu, vv, nx, ny, nz, tmask)

    return ke


# --------------------------------------------------------------------------- #
# Define a function to construct a montly average from a specified
# set of years.

def calc_monthly_ave(fieldname, fieldmask, fieldtype,
                     homedirectory, nemodirectory, filedir, prefix,
                     suffix,
                     start_year=1948, no_years=1, month=1):

    # Import required modules.
    import glob

    # Initialise the output field to zero.
    monthly = 0.0
    counter = 0.0

    # Loop over the number of years and print some shit to the screen
    for k in range(no_years):
        year = start_year + k

        # Find the names of the files.
        possible = ''.join([str(year), '/', filedir, prefix, str(year),
                            str(month).rjust(2, '0'), '??', suffix])
        files = sorted(glob.glob(''.join([homedirectory, nemodirectory,
                                          possible])))

        # Load the current month & year combination into seasonal.
        for m in range(len(files)):
            print(files[m])
            monthly = monthly + load_field(fieldname, '',  '',
                                           files[m], fieldtype)
            counter = counter + 1

    # Divide by the number of files loaded.
    monthly = monthly / counter

    return monthly


# --------------------------------------------------------------------------- #
# Define a function to construct a set of seasonal averages from a specified
# set of years.

def calc_seasonal_ave(fieldname, fieldmask, fieldtype,
                      homedirectory, nemodirectory, filedir, prefix,
                      suffix,
                      start_year=1948, no_years=1, season='JJA'):

    # Import required modules.
    import glob

    # Check that the requested season is one of the actual ones.
    if season not in ('DJF', 'MAM', 'JJA', 'SON', 'JFM', 'AMJ', 'JAS', 'OND'):
        raise Exception('Selected season not in allowed list:',
                        'DJF, MAM, JJA, SON, JFM, AMJ, JAS and OND.')

    # Specify the strings that contain the month to search for
    if season == 'DJF':
        month = ['12??', '01??', '02??']
    elif season == 'MAM':
        month = ['03??', '04??', '05??']
    elif season == 'JJA':
        month = ['06??', '07??', '09??']
    elif season == 'SON':
        month = ['09??', '10??', '11??']
    elif season == 'JFM':
        month = ['01??', '02??', '03??']
    elif season == 'AMJ':
        month = ['04??', '05??', '06??']
    elif season == 'JAS':
        month = ['07??', '08??', '09??']
    elif season == 'OND':
        month = ['10??', '11??', '12??']

    # Initialise the output field to zero.
    seasonal = 0.0
    counter = 0.0

    # Loop over the number of years and print some shit to the screen
    for k in range(no_years):
        for l in range(3):
            if season == 'DJF' and l == 0:
                year = start_year + k - 1
            else:
                year = start_year + k

            # Find the names of the files.
            possible = ''.join([str(year), '/', filedir, prefix, str(year),
                                month[l], suffix])
            files = sorted(glob.glob(''.join([homedirectory, nemodirectory,
                                              possible])))

            # Load the current month & year combination into seasonal.
            for m in range(len(files)):
                print(files[m])
                seasonal = seasonal + load_field(fieldname, '',  '',
                                                 files[m], fieldtype)
                counter = counter + 1

    # Divide by the number of files loaded.
    seasonal = seasonal / counter

    return seasonal


# --------------------------------------------------------------------------- #
# Define a function to construct a set of seasonal averages from a specified
# set of years.

def calc_seasonal_product(field1, field2, fieldmask, fieldtype,
                          homedirectory, nemodirectory, filedir, prefix,
                          suffix,
                          start_year=1948, no_years=1, season='JJA'):

    # Import required modules.
    import glob

    # Check that the requested season is one of the actual ones.
    if season not in ('DJF', 'MAM', 'JJA', 'SON', 'JFM', 'AMJ', 'JAS', 'OND'):
        raise Exception('Selected season not in allowed list:',
                        'DJF, MAM, JJA, SON, JFM, AMJ, JAS and OND.')

    # Specify the strings that contain the month to search for
    if season == 'DJF':
        month = ['12??', '01??', '02??']
    elif season == 'MAM':
        month = ['03??', '04??', '05??']
    elif season == 'JJA':
        month = ['06??', '07??', '09??']
    elif season == 'SON':
        month = ['09??', '10??', '11??']
    elif season == 'JFM':
        month = ['01??', '02??', '03??']
    elif season == 'AMJ':
        month = ['04??', '05??', '06??']
    elif season == 'JAS':
        month = ['07??', '08??', '09??']
    elif season == 'OND':
        month = ['10??', '11??', '12??']

    # Initialise the output field to zero.
    seasonal = 0.0
    counter = 0.0

    # Loop over the number of years and print some shit to the screen
    for k in range(no_years):
        for l in range(3):
            if season == 'DJF' and l == 0:
                year = start_year + k - 1
            else:
                year = start_year + k

            # Find the names of the files.
            possible = ''.join([str(year), '/', filedir, prefix, str(year),
                                month[l], suffix])
            files = sorted(glob.glob(''.join([homedirectory, nemodirectory,
                                              possible])))

            # Load the current month & year combination into seasonal.
            for m in range(len(files)):
                print(files[m])
                seasonal = seasonal + \
                           (load_field(field1, '', '', files[m], fieldtype) *
                            load_field(field2, '', '', files[m], fieldtype))
                counter = counter + 1

    # Divide by the number of files loaded.
    seasonal = seasonal / counter

    return seasonal


# --------------------------------------------------------------------------- #
# Define a function to calculate a scalar product of quantities on U & V points
# that should be on T points.

def calc_uv_on_t(ufield, vfield, nx, ny, nz, tmask):

    # Import required commands.
    from numpy import ndarray
    from numpy.ma import masked_array

    # Define the scalar product field as a masked array full of zeroes.
    scalar_product = masked_array(ndarray(shape=(nx, ny, nz)), mask=1.0-tmask)

    # Calculate the scalar product using 2nd order differences.
    scalar_product.data[0:nx, 0:ny, 0:nz] = 0.5*(ufield[0:nx, 0:ny, 0:nz] +
                                                 ufield[1:nx+1, 0:ny, 0:nz] +
                                                 vfield[0:nx, 0:ny, 0:nz] +
                                                 vfield[0:nx, 1:ny+1, 0:nz])

    return scalar_product


# --------------------------------------------------------------------------- #
# Define a function to calculate the total input of wind power.

def calc_wind_power(utaux, vtauy, taux, tauy, u, v, nx, ny, nz, tmask):

    # Square the velocities.
    utx = u * taux
    vty = v * tauy

    # Calculate the total power input.
    total_power = calc_uv_on_t(utaux, vtauy, nx, ny, nz, tmask)

    # Calculate the power input from the mean wind.
    mean_power = calc_uv_on_t(utx, vty, nx, ny, nz, tmask)

    # Calculate the power input from the wind fluctuations.
    eddy_power = total_power - mean_power

    return total_power, mean_power, eddy_power


# --------------------------------------------------------------------------- #
# Define a function to calculate vorticity from input velocities and grid
# information.

def calc_xi(uvelocity, vvelocity, nx, ny, nz, e1u, e2v, e1f, e2f, fmask):

    # Import required commands.
    from numpy import ndarray
    from numpy.ma import masked_array

    # Define the xi field as a masked_array full of zeroes.
    xi = masked_array(ndarray(shape=(nx+1, ny+1, nz)), mask=1.0-fmask)

    # Calculate xi using 2nd order differences.
    xi.data[1:nx, 1:ny] = (uvelocity[1:nx, 0:ny-1, :] * e1u[1:nx, 0:ny-1, :] +
                           vvelocity[1:nx, 1:ny, :] * e2v[1:nx, 1:ny, :] -
                           uvelocity[1:nx, 1:ny, :] * e1u[1:nx, 1:ny, :] -
                           vvelocity[0:nx-1, 1:ny, :] * e2v[0:nx-1, 1:ny, :]) / \
                          (e1f[1:nx, 1:ny, :] * e2f[1:nx, 1:ny, :])

    # Fill in the grid points at the east and west (ignore north & south as
    # they will be walls).
    xi.data[0:1, 1:ny] = (uvelocity[0:1, 0:ny-1, :] * e1u[0:1, 0:ny-1, :] +
                          vvelocity[0:1, 1:ny, :] * e2v[0:1, 1:ny, :] -
                          uvelocity[0:1, 1:ny, :] * e1u[0:1, 1:ny, :] -
                          vvelocity[-1:, 1:ny, :] * e2v[-1:, 1:ny, :]) / \
                         (e1f[0:1, 1:ny, :] * e2f[0:1, 1:ny, :])
    xi.data[nx:, 1:ny, :] = xi.data[0:1, 1:ny, :]

    return xi


# --------------------------------------------------------------------------- #
# Define a function to construct an annual average from a specified
# set of years.

def extract_field(fieldname, fieldmask, fieldtype,
                  homedirectory, nemodirectory, filedir, prefix,
                  suffix,
                  start_year=1948, no_years=1):

    # Import required modules.
    import glob
    from scipy.io import savemat

    # Initialise the output field to zero.
    counter = 0.0

    # Loop over the number of years and print some shit to the screen
    for k in range(no_years):
        year = start_year + k

        # Find the names of the files.
        possible = ''.join([str(year), '/', filedir, prefix, str(year),
                            '????', suffix])
        files = sorted(glob.glob(''.join([homedirectory, nemodirectory,
                                          possible])))

        # Load the current month & year combination into seasonal.
        for m in range(len(files)):
            print(files[m])
            field = load_field(fieldname, '',  '', files[m], fieldtype)
            counter = counter + 1

            # Save the field to a mat file for matlab-based shenanigans.
            savemat(''.join([files[m][:-9], '_', fieldname, '.mat']),
                    {fieldname: field})

    return counter


# --------------------------------------------------------------------------- #
# Define a function to take the X gradient of a T field.

def grad_dtdx(tfield, e1u, umask, nx ,ny, nz):

    # Import required modules.
    from numpy import ndarray

    # Initialise the dtdx array to have the correct size.
    dtdx = ndarray(shape=(nx+1, ny, nz))

    # Take the X gradient.
    dtdx[1:nx, :, :] = umask[1:nx,:,:] * umask[0:nx-1,:,:] * \
                       (tfield[1:nx, :, :] - tfield[0:nx-1, :, :]) / e1u[1:nx, :, :]

    # Fill in the missing first and last elements.
    dtdx[0, :, :] = umask[0,:,:] * umask[nx-1,:,:] * \
                    (tfield[0, :, :] - tfield[nx-1, :, :]) / e1u[0, :, :]
    dtdx[nx, :, :] = dtdx[0, :, :]

    # Mask the output field.
    dtdx = dtdx * umask

    return dtdx


# --------------------------------------------------------------------------- #
# Define a function to take the Y gradient of a T field.

def grad_dtdy(tfield, e2v, vmask, nx ,ny, nz):

    # Import required modules.
    from numpy import ndarray

    # Initialise the dtdy array to have the correct size.
    dtdy = ndarray(shape=(nx, ny+1, nz))

    # Take the Y gradient.
    dtdy[:, 1:ny, :] = (tfield[:, 1:ny, :] - tfield[:, 0:ny-1, :]) / e2v[:, 1:ny, :]

    # Fill in the missing first and last elements.
    dtdy[:, 0, :] = (tfield[:, 0, :] - tfield[:, ny-1, :]) / e2v[:, 0, :]
    dtdy[:, ny, :] = dtdy[:, 0, :]

    # Mask the output field.
    dtdy = dtdy * vmask

    return dtdy


# --------------------------------------------------------------------------- #
# Define a function to load a field from a netcdf file - eliminates mask (if a
# masked array) and also removes overlaps regions.

def load_field(fieldname, homedirectory, nemodirectory='', filename='',
               fieldtype='', nx=4320, ny=2000, nz=75, bounds=None):

    # Import required commands.
    from numpy import concatenate
    from numpy.ma import is_masked
    from netCDF4 import Dataset

    # Concatenate the locations strings together.
    file_to_load = ''.join([homedirectory, nemodirectory, filename])

    # Open the file for read-only access.
    ncfile = Dataset(file_to_load, mode='r')

    # Load the requested fieldname from the specified file.
    if bounds is None:
        field = ncfile.variables[fieldname][:]
    else:
        field = ncfile.variables[fieldname][bounds]

    # Close the file, for the sake of form.
    ncfile.close()

    # Unmask the fields, because NC4 is doing something fucking stupid.
    if is_masked(field):
        field = field.data

    # Reverse the order of the indices so that its easier to concatenate on
    # extra rows/columns.
    field = field.transpose()

    # Check if the optional argument fieldtype is set, and if it is duplicate
    # the relevant row/column to make the size conform to MITgcm standard.
    if fieldtype == 'T':
        field = field[1:nx+1, 0:]
    elif fieldtype == 'U':
        field = field[1:, 0:]
    elif fieldtype == 'V':
        field = field[1:nx+1, 0:]
        field = concatenate((field[0:, -1:], field[0:, 0:]), axis=1)
    elif fieldtype == 'Z':
        field = field[1:, 0:]
        field = concatenate((field[0:, -1:], field[0:, 0:]), axis=1)

    return field


# --------------------------------------------------------------------------- #
# Define a function to construct a circle in axes coordinates for polar stereographic views.
    
def make_circle():

    # Import required commands - use masked array to preserve masks.
    import numpy as np
    import matplotlib.path as mpath

    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    return circle


# --------------------------------------------------------------------------- #
# Define a function to mask a field that has been previously read in.

def mask_field(field, fieldmask):

    # Import required commands.
    from numpy import ma

    # Define the output masked field with the masked input field and the mask.
    masked_field = ma.masked_array(field*fieldmask, mask=1.0-fieldmask)
    
    # Set the fill_value of the masked_field to be a value that will never appear.
    masked_field.fill_value = -9999.9999

    return masked_field


# --------------------------------------------------------------------------- #
# Define a function to pad a field to east & west for reentrant interpolation.

def pad_field(field, deltalon):

    # Import required commands - use masked array to preserve masks.
    from numpy import ma

    # Pad the field by taking off the longitudinal grid spacing from the first
    # and last longitude.
    padded_field = ma.concatenate((field[0:1]-deltalon,
                                   field,
                                   field[-1:]+deltalon),
                                  axis=0)

    return padded_field


# --------------------------------------------------------------------------- #
# Define a function to pad a field to east & west using reentrant values.

def pad_ocean_field(field):

    # Import required commands - use masked array to preserve masks.
    from numpy import ma

    # Pad the field by adding the last/first zonal values to the beginning/end
    # of the field.
    padded_field = ma.concatenate((field[-1:, :], field, field[0:1, :]),
                                  axis=0)

    return padded_field


# --------------------------------------------------------------------------- #
# Define a function to plot a model field on T points in a south polar stereographic view using cartopy.

def plot_tfield(field, latitude, longitude, minvalue, maxvalue, nocontours, colormapname, dnetxe, gridcolor='white'):

    # Import required commands - use masked array to preserve masks.
    import cartopy
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    # Make circle in axes coordiantes to mask polar stereographic view.
    circle = make_circle()

    # Make a landmask that is 1 in land and zero in ocean.
    #landmask = np.ma.masked_where(field.mask ==0.,field.mask) 

    # Create a new figure window.
    plt.figure(figsize=(12.0, 12.0))

    # Create the axis with the Greenwich meridian point northwest.
    axisname = plt.axes(projection=cartopy.crs.SouthPolarStereo(central_longitude=135))

    # Geoaxes ignoes facecolor, so make the background black this way.
    axisname.background_patch.set_facecolor('k')

    # Try and use set_bad to make land black.
    #blarp = plt.get_cmap(colormapname)
    #blarp.set_bad(color='black', alpha=1.0)

    # Set the boundary to be the defined circle.
    axisname.set_boundary(circle, transform=axisname.transAxes, use_as_clip_path=True)

    # Draw the contour plot of field.
    cs = axisname.contourf(longitude, latitude, 
                           field, levels=np.linspace(minvalue, maxvalue, nocontours),
                           transform=cartopy.crs.PlateCarree(),
                           cmap=colormapname, vmin=minvalue, vmax=maxvalue,
                           extend=dnetxe)

    # Use the fieldmask to draw solid black where land is.
    #cs_mask = axisname.contourf(longitude, latitude, landmask, [1.0, 1.01],
    #                            colors='black',
    #                            transform=cartopy.crs.PlateCarree())

    # Put the northern boundary at -35oS.
    axisname.set_extent([-45., 315., -90., -35.], cartopy.crs.PlateCarree())

    # Add gridlines.
    gl = axisname.gridlines(crs=cartopy.crs.PlateCarree(), linewidth=1, color=gridcolor, alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator([-45., 0., 45., 90., 135., 180., 225., 270, 315.])
    gl.ylocator = mticker.FixedLocator([-75., -60., -45., -30.])
    gl.n_steps = 360

    # Add a colour bar.
    colorbarname = plt.colorbar(cs, ticks=np.linspace(minvalue, maxvalue, 5), fraction=0.046, pad=0.04)
    colorbarname.ax.tick_params(labelsize=18)
    colorbarname.outline.set_linewidth(1)

    return axisname, colorbarname


# --------------------------------------------------------------------------- #
# Define a function to calculate the heat transort and provide Reynolds averaged components.
    
def reynolds_average_vt(vt, v, t, e1v, e3t, e3v, rho0= 1026.0, cp=3991.86795711963):

    # Import required commands - use masked array to preserve masks.
    from numpy import sum
    
    # Vertically and zonally integrate to get the total heat transport.
    total = rho0 * cp * sum(e1v[:, :, None] * vt * e3v, axis=(0, 2))

    # Use the mean V & T to calculate the mean heat transport : mean = zonal + standing.
    mean = rho0 * cp * sum(e1v[:, :, None] * v * t * e3v, axis=(0, 2))

    # Calculate the zonal mean T.
    zonalt = sum(e1v[:, :, None] * e3v * t, axis=0) / sum(e1v[:, :, None] * e3v, axis=0)

    # Calculate the zonal mean V.
    zonalv = sum(e1v[:, :, None] * e3v * v, axis=0) / sum(e1v[:, :, None] * e3v, axis=0)

    # Use the zonal mean V & T to calculate the heat transport by the zonal mean.
    zonal = rho0 * cp * sum( e1v[:, :, None] * zonalv[None, :, :] * zonalt[None, :, :] * e3v, axis=(0, 2))

    # Calculate the standing eddy pattern.
    vtstar = (v - zonalv[None, :, :]) * (t - zonalt[None, :, :])

    # Use the standing eddy component to calculate the standing eddy heat transport.
    standing = rho0 * cp * sum(e1v[:, :, None] * vtstar * e3v, axis=(0, 2))

    # Calculate the covariance of V & T.
    vtprime = vt - v*t

    # Use the covariance of V & T to calculate the transient heat transport.
    transient = rho0 * cp * sum(e1v[:, :, None] * vtprime * e3v, axis=(0, 2))

    return {'total' : total, 'mean' : mean, 'zonal' : zonal, 'standing' : standing, 'transient' : transient}
    

# --------------------------------------------------------------------------- #
