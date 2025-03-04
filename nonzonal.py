import numpy as np


# Define filename and path to save data
filename = ... # TODO: specify filename (string)
path2data = ... # TODO: change this to the path to the data (string)

# Compute bottom friction coefficient
mu_days = ... # TODO: specify mu^{-1} in days (float)
if mu_days == 0:
    mu = 0
else:
    mu = 1/(mu_days*86400)

# Angle of velocity vector
angle_UV = ... # TODO: specify angle of velocity vector in degrees (float)

g = 9.81



def c_im_max(kappa,theta,alpha_x,alpha_y,mu,f0,beta,U,V,rho1,rho2,H1,H2):
    """
    Computes the growth rate of the most unstable mode for the 2-layer QG model

    kappa, theta: wavenumber magnitude and angle
    alpha_x, alpha_y: x and y components of bottom slope
    mu: inverse frictional timescale
    f0: Coriolis parameter
    beta: meridional gradient of Coriolis parameter
    U, V: background velocity shear
    rho1, rho2: density of the two layers
    H1, H2: depth of the two layers
    """
    # Calculating relevant parameters
    gp = g*(rho2-rho1)/rho2
    F1 = f0**2/(gp*H1)
    F2 = f0**2/(gp*H2)

    # Coordinate transformation: unit vectors in tangential and normal direction of wavenumber vector
    cos_theta = np.around(np.cos(theta),decimals=15) # to avoid roundoff errors for sin(180), cos(90) and cos(270)
    sin_theta = np.around(np.sin(theta),decimals=15)
    U_new = U*cos_theta + V*sin_theta
    beta_new = beta*cos_theta
    S = f0/H2 * (alpha_y*cos_theta - alpha_x*sin_theta)
    
    # Solving the equation
    a1 = kappa**2 * (kappa**2 + F1 + F2)
    a2 = - U_new*kappa**2 * (kappa**2 + 2*F2) + beta_new * (2*kappa**2 + F1 + F2) + (kappa**2 + F1) * (S + 1j*mu*kappa)
    a3 = (-U_new*kappa**2 + beta_new) * (-F2*U_new + beta_new + S + 1j*mu*kappa)
    D = a2**2 - 4*a1*a3

    sol1 = (-a2 + np.sqrt(D))/(2*a1)
    sol2 = (-a2 - np.sqrt(D))/(2*a1)

    # Computing the imaginary parts of the solutions
    cim1 = np.imag(sol1)
    cim2 = np.imag(sol2)
    cre1 = np.real(sol1)
    cre2 = np.real(sol2)

    # If growth rates are negative, make them NaN (no growth)
    if cim1<=0:
        cim1 = np.nan
    if cim2<=0:
        cim2 = np.nan

    max_growth = np.nanmax((cim1, cim2))
    if np.isnan(max_growth):
        return np.nan, np.nan
    elif max_growth == cim1:
        return kappa*cim1, cre1
    else:
        return kappa*cim2, cre2 # return max growth rate & corresponding phase speed



def mostUnstable(kappa_range, alpha_x, alpha_y, mu, f0, beta, U, V, rho1, rho2, H1, H2, F):
    """
    Returns the instability properties for a given range of wavenumbers
    - most unstable wavenumber normalized by deformation radius
    - maximum growth rate in 1/days
    - phase speed in cm/s
    - angle of phase velocity vector in degrees
    """
    theta_range = np.deg2rad(np.arange(0,361,2))
    c_range = np.array([[c_im_max(kappa, theta, alpha_x, alpha_y, mu, f0, beta, U, V, rho1, rho2, H1, H2) 
                            for theta in theta_range] for kappa in kappa_range])
    growth_rate = c_range[:,:,0]  # growth rate
    phase_speed = c_range[:,:,1]  # phase speed

    if np.isnan(growth_rate).all():
        return np.nan, np.nan, np.nan, np.nan
    else:
        ind = np.unravel_index(np.nanargmax(growth_rate, axis=None), growth_rate.shape)
        if phase_speed[ind] >= 0:
            c_ang = 0
        elif theta_range[ind[1]] < np.pi:
            c_ang = np.pi
        else:
            c_ang = -np.pi
        return kappa_range[ind[0]]/np.sqrt(F), growth_rate[ind]*86400, phase_speed[ind]*1e2, np.rad2deg(theta_range[ind[1]]+c_ang)
    


# Define parameters
H1 = 1000
H2 = 4000
rho1 = 1027.5
rho2 = 1028
f0 = 1e-4
beta = 1e-11

gp = g*(rho2-rho1)/rho2
F1 = f0**2/(gp*H1)
F2 = f0**2/(gp*H2)
F = F1 + F2

U = 0.04*np.cos(angle_UV*np.pi/180)
V = 0.04*np.sin(angle_UV*np.pi/180)



# For a variety of bottom slope vectors (different magnitudes/orientations), find the most unstable mode
k_range = np.arange(0.2,1.31,0.01)*np.sqrt(F)
alpha_range = np.arange(0,3.1e-3,5e-5)
phi_range = np.deg2rad(np.arange(0,361,1))
k_inj = np.zeros((len(alpha_range),len(phi_range)))
growth_rate_inj = np.zeros((len(alpha_range),len(phi_range)))
phase_speed_inj = np.zeros((len(alpha_range),len(phi_range)))
theta_c = np.zeros((len(alpha_range),len(phi_range)))

for i in range(len(alpha_range)):
    for j in range(len(phi_range)):
        alpha_x = alpha_range[i]*np.cos(phi_range[j])
        alpha_y = alpha_range[i]*np.sin(phi_range[j])
        k_inj[i,j],growth_rate_inj[i,j],phase_speed_inj[i,j],theta_c[i,j] = mostUnstable(k_range, alpha_x, alpha_y, mu, f0, beta, U, V, rho1, rho2, H1, H2, F)


np.save(path2data+filename+'_k.npy', k_inj)
np.save(path2data+filename+'_sigma.npy', growth_rate_inj)
np.save(path2data+filename+'_cr.npy', phase_speed_inj)
np.save(path2data+filename+'_theta_cr.npy', theta_c)


print("Run with 1/mu = %i days, angle_UV = %i degrees finished"%(mu_days,angle_UV))
