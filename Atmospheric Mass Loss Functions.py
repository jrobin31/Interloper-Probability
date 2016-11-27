# The purpose of this code is to break the components of the atmospheric mass loss calculation into functions
import numpy as np
import matplotlib.pyplot as pp
from astropy.io import ascii
from astropy.table import Table




# Global parameters
_bol_luminosity = 1e46

# The first function we will define is the emission spectrum of the AGN
def emission_spec():
    richards = ascii.read('/Users/Jenna/Documents/Matt Mechtley Quasar Research/Richards_Mean_Quasar_Data.txt')
    wavelength = 3e14 / 10**richards['LogF']
    flux = richards['All']
    flux = 10**flux
    flux = flux / (10**richards['LogF'])
    my_table = Table([wavelength, flux], names=['wavelength', 'flux']) 
    my_table['wavelength'].unit = 'microns'
    my_table['flux'].unit = 'erg / s / Hz'
    print(my_table)
    return my_table
    # Need wavelength and flux columns --> wavelength v. flux before dust & wavelength v. flux after dust

# Calculate flux at planet using luminosity
def spec_in_flux(spec_in_luminosity, distance):
    spec_in_luminosity['flux'] = spec_in_luminosity['flux'] / (4*np.pi*distance**2)
    return spec_in_luminosity

# The next function we will define is the luminosity (the normalization of the brightness of the spectrum)
def luminosity():
    # Essentially a constant value that we can plug in
    # Multiply table column by constant
    # Bolometric correction at 5100 angstroms = 12.17

    return

# The third function we will define is the distance from the quasar as a function of the orbital parameters
def r_phi_z_position(current_position, current_velocity, ang_momentum, q, delta_time):
    # Parameters will include eccentricity of orbit and angular momentum
    radius, phi, height = current_position # Vector unpacking
    phi_dot = ang_momentum / radius ** 2
    vel_r_dot = ang_momentum ** 2 / radius ** 3 \
               - radius / (radius ** 2 + height ** 2)
    vel_z_dot = -height / q / (radius ** 2 + height ** 2)

    # Calculate all new positions from old velocity, then calculate new velocity
    r_phi_z_dot = np.array((current_velocity[0], phi_dot, current_velocity[1] / q))
    new_position = current_position + r_phi_z_dot * delta_time
    # Calculates the change in the position with the time step
    new_velocity = current_velocity + np.array((vel_r_dot, vel_z_dot)) * delta_time

    return new_position, new_velocity

# The fourth function we will define is the metallicity of the host star system
def metallicity_host():
    # Z component of the original calculation - do we need to manipulate the original equation to get this?
    # We need to talk to Patrick Young on this one
    # This will consider photochemistry as well - talk to Ariel Anbar 

    return

def torus_obscuration(position, luminosity, r_in=20, r_out=600, tau_e=10, beta=1.1, gamma=1):
    # Units for radius: AU; min. 20, max. 600
    theta_h = np.arccos((1 + 3*luminosity / 10**42.65)**(1 - 2*0.44))
    # Values for theta_h come from Lusso et al paper and generate the dashed line from figure 19 right panel
    equatorial_density = r_out**(1-beta) / (1-beta) - r_in**(1-beta) / (1-beta)
    density_constant = tau_e / equatorial_density
    x_y_z_position = np.array((position[0]*np.cos(position[1]), position[0]*np.sin(position[1]), position[2]))
    polar_angle = np.dot(x_y_z_position, (0,0,1)) / np.sqrt(np.dot(x_y_z_position, x_y_z_position))

    if np.arccos(np.abs(polar_angle)) < theta_h:
        return 0
    else:
        return 0.61*density_constant*equatorial_density*np.exp(-gamma*np.abs(polar_angle))

# The last function we will define is the dust obscurity as a function of where the planet lies in the taurus
def dust_obscurity(spectrum, a_v):
    # dust_obscurity will essentially calculate the angle between the planet and the axis of the quasar
    # k_lam is only valid for UV wavelengths (0.12 - 0.63 microns)
    spectrum = spectrum.copy()
    lam = spectrum['wavelength']
    k_lam = np.zeros_like(lam)
    uv_opt = (lam < 0.63) & (lam > 0.12)
    k_lam[uv_opt] = 2.569*(-2.156+(1.509/lam[uv_opt])-(0.198/lam[uv_opt]**2)+(0.011/lam[uv_opt]**3))+4.05
    opt_nir = (lam >= 0.63) & (lam < 2.2)
    k_lam[opt_nir] = 2.659*(-1.857+(1.040/lam[opt_nir]))+4.05
    # X-ray UV eventually
    # 4.05 is the Calzetti R_V value
    a_lam = k_lam * a_v / 4.05
    # Converting attentuation in magnitudes to attenuation in flux
    spectrum['flux'] = spectrum['flux'] * 10**(-0.4 * a_lam)
    return spectrum

def main_loop():
    step_size = 1e-3
    steps = int(10 / step_size)
    # Initial conditions
    r_phi_z_0 = np.array((0.31, 0.0, 0.17))
    vel_0 = np.array((0.0, 0.0))
    q = 0.9
    # Angular momentum of the orbit
    Lz = 0.2

    spectrum = emission_spec()

    # Calculate velocities at t=1/2 (leapfrog integration)
    pos_half, vel_half = r_phi_z_position(r_phi_z_0, vel_0, Lz, q, step_size / 2.0)
    # 2.0 rather than 2 because we want a floating point value rather than an integer

    current_position = r_phi_z_0
    current_velocity = vel_half

    all_positions = np.zeros((steps, 3))
    all_avs = np.zeros(steps)

    for step in range(steps):
        all_positions[step, :] = current_position
        current_position, current_velocity = r_phi_z_position(current_position, current_velocity, Lz, q, step_size)
        all_avs[step] = torus_obscuration(current_position, _bol_luminosity)
        red_spec = dust_obscurity(spectrum, all_avs[step])
        flux_spec = spec_in_flux(red_spec, np.sqrt(current_position[0]**2 + current_position[2]**2)) # Need to put this into real units
        # Flux_spec is the thing that actually interacts with the atmosphere

    pp.plot(all_positions[:, 0], all_positions[:, 2])
    pp.xlabel('Radius')
    pp.ylabel('Disk Height')
    # No title on plots if you're putting them in a journal
    pp.show()

# Plot y vs. x
    x = all_positions[:, 0]*np.cos(all_positions[:, 1])
    y = all_positions[:, 0]*np.sin(all_positions[:, 1])
    pp.plot(x, y, label='$\~L_z = %0.2g$, $q = %0.2g$' % (Lz, q))
    pp.xlabel('$x$')
    pp.ylabel('$y$')
    pp.legend(loc='upper right')
    pp.title('X/Y Position Plot')
    pp.show()

# Plot av vs. time
    pp.plot(range(steps), all_avs)
    pp.xlabel('Time step')
    pp.ylabel('$A_V$')
    pp.show()

# Find zero crossings (in terms of radius and radial velocity)
# Zero crossings are where the star crosses the zero point in the z-plane (when it moves down); tells us if the star is
# moving toward or away from the center
#up_zc_mask = np.diff(np.sign(heights), 1) > 0

# START EDITING HERE
#plot_title = r'$\~v_{r,\tau=0} = %0.2g$, $\~v_{z,\tau=0} = %0.2g$,' \
             #r'$\~r_{\tau=0} = %0.2g$, $\~z_{\tau=0} = %0.2g$' % \
             #(vel_r0, vel_z0, rad0, height0) #r_phi_z_0, vel_half

# Plot upward zero crossings, r vs. v_r
#ax.plot(radii[1:][up_zc_mask], vel_rs[1:][up_zc_mask], 'o',
        #label='$\~L_z = %0.2g$, $q = %0.2g$' % (Lz, q)) # Might need to fix this too... Check once other calculations
# are complete
#ax.set_xlabel('$\~r$')
#ax.set_ylabel('$\~v_r$')
#ax.legend(loc='lower right')
#ax.set_title(plot_title)
#pp.show()


if __name__ == "__main__":
    main_loop()
    # spectrum = emission_spec()
    # pp.loglog(spectrum['wavelength'], spectrum['flux'], color="blue")
    # pp.xlabel('Wavelength ($\mu$m)')
    # pp.ylabel('Flux (erg / s / Hz)')
    # spectrum = dust_obscurity(spectrum, 2.0) # a_v here is a free variable
    # pp.loglog(spectrum['wavelength'], spectrum['flux'], color="red")
    # pp.xlabel('Wavelength ($\mu$m)')
    # pp.ylabel('Flux (erg / s / Hz)')
    # spectrum = spec_in_flux(spectrum, distance=10) # distance here is also a free variable
    # pp.loglog(spectrum['wavelength'], spectrum['flux'], color="green")
    # pp.xlabel('Wavelength ($\mu$m)')
    # pp.ylabel(r'Flux ($f_\nu$, arbitrary)')
    # pp.show() 