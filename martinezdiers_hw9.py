#Nina Martinez Diers
#Adapted from MD simulation provided by Vianney Gimenez-Pinto with xyz trajecotry file output, 
#which was an example from pythoninclemistry.org

#2-D with lj_force, accel, velocity, position functions updated
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann
mass_of_argon = 39.948 # amu

def lj_force(r, epsilon, sigma): #unchanged from original
    """
    Implementation of the Lennard-Jones potential 
    to calculate the force of the interaction.
    
    Parameters
    ----------
    r: float
        Distance between two particles (Å)
    epsilon: float 
        Potential energy at the equilibrium bond 
        length (eV)
    sigma: float 
        Distance at which the potential energy is 
        zero (Å)
    
    Returns
    -------
    float
        Force of the van der Waals interaction (eV/Å)
    """
    return 48 * epsilon * np.power(
        sigma, 12) / np.power(
        r, 13) - 24 * epsilon * np.power(
        sigma, 6) / np.power(r, 7)


def init_velocity(T, number_of_particles):
    """
    Initialise the velocities for a series of 
    particles.
    
    Parameters
    ----------
    T: float
        Temperature of the system at 
        initialisation (K)
    number_of_particles: int
        Number of particles in the system
    
    Returns
    -------
    ndarray of floats
        Initial x-components of velocities for a series of 
        particles (eVs/Åamu)

    ndarray of floats
        Initial y-components of velocities for a series of 
        particles (eVs/Åamu)
    """
    Rx = np.random.rand(number_of_particles) - 0.5
    #Ry = np.zeros(number_of_particles) #for 1-D case
    Ry = np.random.rand(number_of_particles) - 0.5 #for 2-D case
    return Rx * np.sqrt(Boltzmann * T / (
        mass_of_argon * 1.602e-19)), Ry * np.sqrt(Boltzmann * T / (
        mass_of_argon * 1.602e-19))


def get_accelerations(xpositions,ypositions): #takes 2 arguments instead of 1 to calculate both components of acceleration
    """
    Calculate the acceleration on each particle
    as a  result of each other particle. 
    N.B. We use the Python convention of 
    numbering from 0.
    
    Parameters
    ----------
    xpositions: ndarray of floats
        The positions, in the x dimension, 
        for all of the particles
    ypositions: ndarray of floats
        The positions, in the y dimension, 
        for all of the particles
        
    Returns
    -------
    ndarray of floats
        The x-component of acceleration on each
        particle (eV/Åamu)
    ndarray of floats
        The y-component of acceleration on each
        particle (eV/Åamu)
    """
    accel_x = np.zeros((xpositions.size, xpositions.size))
    accel_y = np.zeros((ypositions.size, ypositions.size))
    for i in range(0, xpositions.size - 1):
        for j in range(i + 1, xpositions.size):
            r_x = xpositions[j] - xpositions[i]
            r_y = ypositions[j] - ypositions[i]
            rmag = np.sqrt((r_x * r_x) + (r_y * r_y))
            force_scalar = lj_force(rmag, 0.0103, 3.4)
            force_x = force_scalar * r_x / rmag
            accel_x[i, j] = force_x / mass_of_argon
            accel_x[j, i] = - force_x / mass_of_argon
            
            force_y = force_scalar * r_y / rmag
            accel_y[i, j] = force_y / mass_of_argon
            accel_y[j, i] = - force_y / mass_of_argon
    return np.sum(accel_x, axis=0), np.sum(accel_y, axis=0)

def update_pos(x, y, vx, vy, ax, ay, dt):
    """
    Update the particle positions in two dimensions.
    
    Parameters
    ----------
    x: ndarray of floats
        The x-positions of the particles in 
        two dimensions
    y: ndarray of floats
        The y-positions of the particles in 
        two dimensions
    vx: ndarray of floats
        The x-component of velocities of the 
        particles in two dimensions
    vy: ndarray of floats
        The y-component of velocities of the 
        particles in two dimensions
    ax: ndarray of floats
        The x-component of accelerations of the 
        particles in two dimensions.
    ay: ndarray of floats
        The y-component of accelerations of the 
        particles in two dimensions.
    dt: float
        The timestep length
    
    Returns
    -------
    ndarray of floats:
        New x-components of positions of the 
        particles in two dimensions
    ndarray of floats:
        New y-components of positions of the 
        particles in two dimensions
    """
    return x + vx * dt + 0.5 * ax * dt * dt, y + vy * dt + 0.5 * ay * dt * dt

def update_velo(vx, vy, ax, ay, ax1, ay1, dt):
    """
    Update the particle velocities in two dimensions.
    
    Parameters
    ----------
    vx: ndarray of floats
        The x-component of velocities of the 
        particles in two dimensions
    vy: ndarray of floats
        The y-component of velocities of the 
        particles in two dimensions
    ax: ndarray of floats
        The x-component of accelerations of the 
        particles in two dimensionsat the previous 
        timestep (eV/Åamu).
    ay: ndarray of floats
        The y-component of accelerations of the 
        particles in two dimensionsat the previous 
        timestep (eV/Åamu).
    a1x: ndarray of floats
        The x-component of accelerations of the 
        particles in two dimensions at the current 
        timestep (eV/Åamu)
    a1y: ndarray of floats
        The y-component of accelerations of the 
        particles in two dimensions at the current 
        timestep (eV/Åamu)
    dt: float
        The timestep length
    
    Returns
    -------
    ndarray of floats:
        New x-component of velocities of the particles 
        in two dimensions (eVs/Åamu)
    ndarray of floats:
        New y-component of velocities of the particles 
        in two dimensions (eVs/Åamu)
    """
    return vx + 0.5 * (ax + ax1) * dt, vy + 0.5 * (ay + ay1) * dt

def run_md(dt, number_of_steps, snaptime, initial_temp, x, y, number_of_atoms):
    """
    Run a MD simulation in two dimensions.
    
    Parameters
    ----------
    dt: float
        The timestep length (s)
    number_of_steps: int
        Number of iterations in the simulation
    snaptime: int
        Sampling for MD simulation frames
    initial_temp: float
        Temperature of the system at 
        initialisation (K)
    x: ndarray of floats
        Initial x-positions of the particles in 
        two dimensions
    y: ndarray of floats
        Initial y-positions of the particles in 
        two dimensions
    number_of_atoms: int
        Number of atoms in the simulation
        
    Returns
    -------
    ndarray of floats:
        The x-component of positions for all of the particles 
        throughout the simulation (Å)
    ndarray of floats:
        The x-component of positions for all of the particles 
        throughout the simulation (Å)
    """

    file = open('traj.xyz', 'w')
    xpositions = np.zeros((number_of_steps, number_of_atoms))
    ypositions = np.zeros((number_of_steps, number_of_atoms))
    vx,vy = init_velocity(initial_temp, number_of_atoms)
    ax,ay = get_accelerations(x,y)
    for i in range(number_of_steps):
        x,y = update_pos(x,y, vx,vy, ax,ay, dt)
        ax1,ay1 = get_accelerations(x,y)
        vx,vy = update_velo(vx,vy, ax,ay, ax1,ay1, dt)
        ax = np.array(ax1)
        ay = np.array(ay1)
        xpositions[i, :] = x
        ypositions[i, :] = y
        if i%snaptime == 0:
            file.write(str(number_of_atoms)+"\n")
            file.write("#\n")
            for atom in range(number_of_atoms):
               # print(atom)
                file.write("A\t" + str(x[atom])+"\t"+ str(y[atom])+ "\t0.0\n" )
    return xpositions,ypositions



number_of_atoms=10
x=1+5*np.arange(number_of_atoms)
y=1+5*np.arange(number_of_atoms)
sim_posx, sim_posy = run_md(0.1, 10000,100, 300, x, y, number_of_atoms)

#plot trajectories of particles
%matplotlib inline
for i in range(sim_posx.shape[1]):
    plt.plot(sim_posx[:, i], sim_posy[:,i], '.', label='atom {}'.format(i))
plt.xlabel(r'$x$-Position/Å')
plt.ylabel(r'$y$-Position/Å')
plt.legend(frameon=False)
plt.show()