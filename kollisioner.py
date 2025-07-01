"""
Elastic collisions.
Energy transfer from light projectile to heavy target, using equation 2.6b 
from "Backscattering Spectrometry" by Chu, Mayer and Nicolet.
Program written in June 2020, November 2020, updated plotting December 2021
Updated for mechanics course February 2025
"""
%matplotlib qt
import numpy as np
import matplotlib.pyplot as plt


# Decide if you want to have one common or three separate figures
common = False


# Part I
# Use the equation mentioned above


# Range for the mass ratio M1/M2, M1 = projectile, M2 = target
# Light projectiles on heavy targets:
mass_ratios = np.linspace(0, 1, 51)
# Range for the scattering angle of the projectile
Theta_values = np.linspace(0, np.pi, 91)
# Create a grid with mass ratios and angles
X, Y = np.meshgrid(mass_ratios, Theta_values)


# Energy transfer: 1 - Ratio between the projectile energy after and before the collision
Z = 1 - (( X * np.cos(Y) + np.sqrt(1 - X**2 * np.sin(Y)**2)) / (1+X))**2


# Part II
# Simulating particles with different impact parameters for the collision

# Initial velocity in x-direction.
v0 = np.matrix([1,0])
# The individual radii don't matter, only the ratio between their sum and 
# the impact parameter.
R1 = 0.5
R2 = 0.5


# Simulate lots of projectiles.
number_of_simulations = 50000
data = np.zeros([number_of_simulations, 4])
n = 0
while n < number_of_simulations:
    # Random start position in the xy-plane
    random_pair = np.random.uniform(-(R1+R2), (R1+R2), 2)
    # Impact parameter.
    d = np.sqrt((random_pair**2).sum())
    if d > (R1+R2):
        # The projectile missed the target
        continue
    alpha = np.arcsin(d/(R1+R2))   # Always > 0
    # v0 in a radial/tangential coordinate system just before the collision.
    A = np.matrix([[np.cos(alpha), np.sin(alpha)], [-np.sin(alpha), np.cos(alpha)]])
    v0_rt = v0 * A
    # A random mass ratio from mass_ratios.
    x = mass_ratios.min() + np.random.rand()*(mass_ratios.max()-mass_ratios.min())
    # Central collision for the radial component.
    B = np.matrix([[(x-1)/(x+1), 0], [0, 1]])
    v1_rt = B * v0_rt.transpose()
    # v1 = (v1_rt.transpose() * A.transpose()).transpose()
    v1 = A * v1_rt
    Theta = np.arccos(v1[0]/np.sqrt((np.array(v1)**2).sum()))[0,0]
    z = 1 - (np.array(v1_rt)**2).sum()
    data[n] = [x, Theta, z, d]
    n = n+1



if common:
    # Plotting, all in one figure
    # Analytic solution
    fig = plt.figure(figsize=(18,10))
    ax = fig.add_subplot(221, projection='3d')
    ax.set_title('Energy transfer from projectile to target, analytic')
    ax.set_xlabel('Mass ratio M(projectile)/M(target)')
    ax.set_ylabel('Scattering angle')
    ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
    ax.set_zlabel('Relative energy transfer')
    ax.plot_surface(X, Y*180/np.pi, Z, cmap='rainbow',  rcount=len(mass_ratios), ccount=len(Theta_values))
    
    # Simulation, scattering angle vs impact parameter
    ax = fig.add_subplot(223, projection='3d')
    ax.set_title('Scattering angle vs impact parameter, simulation')
    ax.set_xlabel('Mass ratio M(projectile)/M(target)')
    ax.set_ylabel('Impact parameter')
    #ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
    ax.set_zlabel('Scattering angle')
    ax.scatter(data[:,0], data[:,3], data[:,1]*180/np.pi, s=10, c=data[:,1], cmap='rainbow')
    
    # Simulation, projectile energy vs scattering angle 
    ax = fig.add_subplot(224, projection='3d')
    ax.set_title('Energy transfer from projectile to target, simulation')
    ax.set_xlabel('Mass ratio M(projectile)/M(target)')
    ax.set_ylabel('Scattering angle')
    ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
    ax.set_zlabel('Relative energy transfer')
    ax.scatter(data[:,0], data[:,1]*180/np.pi, data[:,2], s=10, c=data[:,2], cmap='rainbow')
    fig.show()

else:
    # Plotting, three figures
    # Analytic solution
    fig = plt.figure(figsize=(18,10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title('Energy transfer from projectile to target, analytic')
    ax.set_xlabel('Mass ratio M(projectile)/M(target)')
    ax.set_ylabel('Scattering angle')
    ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
    ax.set_zlabel('Relative energy transfer')
    ax.plot_surface(X, Y*180/np.pi, Z, cmap='rainbow',  rcount=len(mass_ratios), ccount=len(Theta_values))
    fig.show()
    
    # Simulation, scattering angle vs impact parameter
    fig = plt.figure(figsize=(18,10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title('Scattering angle vs impact parameter, simulation')
    ax.set_xlabel('Mass ratio M(projectile)/M(target)')
    ax.set_ylabel('Impact parameter')
    #ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
    ax.set_zlabel('Scattering angle')
    ax.scatter(data[:,0], data[:,3], data[:,1]*180/np.pi, s=10, c=data[:,1], cmap='rainbow')
    
    # Simulation, projectile energy vs scattering angle 
    fig = plt.figure(figsize=(18,10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title('Energy transfer from projectile to target, simulation')
    ax.set_xlabel('Mass ratio M(projectile)/M(target)')
    ax.set_ylabel('Scattering angle')
    ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
    ax.set_zlabel('Relative energy transfer')
    ax.scatter(data[:,0], data[:,1]*180/np.pi, data[:,2], s=10, c=data[:,2], cmap='rainbow')
    fig.show()



# The following is some analysis of the probability for the different scattering angles,
# based on the simulation above. Uncomment it if you want to test...

# # Probabilities for different scattering angles, part 1
# fig = plt.figure(figsize=(12,6))
# plt.title('Scattering probabilities (3d). Light projectiles, mass ratio < 0.05')
# plt.hist(data[data[:,0]<0.05][:,1]*180/np.pi, bins=[n * 15 for n in range(13)])
# plt.xticks([0, 30, 60, 90, 120, 150, 180])
# plt.xlabel('Scattering angle')
# plt.ylabel('Counts')
# plt.show()

# # Probabilities for different scattering angles, part 2
# fig = plt.figure(figsize=(12,6))
# plt.title('Scattering probabilities (3d). Heavy projectiles, mass ratio > 0.98')
# plt.hist(data[data[:,0]>0.98][:,1]*180/np.pi, bins=[n * 15 for n in range(13)])
# plt.xticks([0, 30, 60, 90, 120, 150, 180])
# plt.xlabel('Scattering angle')
# plt.ylabel('Counts')
# plt.show()


"""
Toves histogram
"""
X, Y = np.meshgrid(mass_ratios, Theta_values)
Z = 1 - (( X * np.cos(Y) + np.sqrt(1 - X**2 * np.sin(Y)**2)) / (1+X))**2

hist_values, _, _ = np.histogram2d(data[:,0], data[:,1]*180/np.pi, bins=[X.shape[0], Y.shape[1]])
hist_norm = hist_values / np.max(hist_values)  # Normalize between 0 and 1

fig = plt.figure(figsize=(18,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Energy transfer from projectile to target, analytic')
ax.set_xlabel('Mass ratio M(projectile)/M(target)')
ax.set_ylabel('Scattering angle')
ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
ax.set_zlabel('Relative energy transfer')
ax.plot_surface(X, Y*180/np.pi, Z, facecolors=plt.cm.viridis(hist_norm), edgecolor='k')

mappable = plt.cm.ScalarMappable(cmap="viridis")
mappable.set_array([])  
fig.colorbar(mappable, ax=ax, label="Histogram Frequency")

fig.show()