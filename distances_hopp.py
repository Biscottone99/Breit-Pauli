import numpy as np
import matplotlib.pyplot as plt

l = 5
c = 2
alpha = np.linspace(0.01, np.pi/2, 100)  # eviti divisione per zero
distance = np.zeros((len(alpha), 2))
hopping = np.zeros((len(alpha), 2))
norm_cdot = np.zeros((len(alpha), 3))
for idx, a in enumerate(alpha):
    r = np.sqrt((l**2 - c**2) / (2 * (1 - np.cos(a))))
    site1 = [r * np.cos(np.pi/2), r * np.sin(np.pi/2), 0]
    site2 = [r * np.cos(np.pi/2 - a), r * np.sin(np.pi/2 - a), c]
    site3 = [r * np.cos(np.pi/2 - 2*a), r * np.sin(np.pi/2 - 2*a), 2*c]
    site4 = [r * np.cos(np.pi/2 - 3*a), r * np.sin(np.pi/2 - 3*a), 3*c]

    d13 = np.sqrt((site1[0] - site3[0])**2 + (site1[1] - site3[1])**2 + (site1[2] - site3[2])**2)
    d14 = np.sqrt((site1[0] - site4[0])**2 + (site1[1] - site4[1])**2 + (site1[2] - site4[2])**2)

    hop13 = 0.1 * np.exp(-d13)
    hop14 = 0.1 * np.exp(-d14)
    cdot13= np.cross(site4, site2)
    cdot14= np.cross(site3, site4)
    cdot12= np.cross(site4, site1)
    norm_cdot13 = np.linalg.norm(cdot13)
    norm_cdot14 = np.linalg.norm(cdot14)
    norm_cdot12 = np.linalg.norm(cdot12)
    norm_cdot[idx, 0] = norm_cdot13
    norm_cdot[idx, 1] = norm_cdot14
    norm_cdot[idx, 2] = norm_cdot12
    distance[idx, 0] = d13
    distance[idx, 1] = d14
    hopping[idx, 0] = hop13
    hopping[idx, 1] = hop14

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)

# --- Distanze ---
axes[0].plot(alpha*180/np.pi, distance[:,0], label='Distance 1-3')
axes[0].plot(alpha*180/np.pi, distance[:,1], label='Distance 1-4')
axes[0].set_xlabel(r'$\alpha$ (degrees)')
axes[0].set_ylabel('Distance')
axes[0].set_title('Distances vs Alpha')
axes[0].set_xlim(0, 90)
axes[0].legend(frameon=False)

# --- Hopping ---
axes[1].plot(alpha*180/np.pi, hopping[:,0], label='Hopping 1-3')
axes[1].plot(alpha*180/np.pi, hopping[:,1], label='Hopping 1-4')
axes[1].set_xlabel(r'$\alpha$ (degrees)')
axes[1].set_ylabel('Hopping')
axes[1].set_title('Hopping vs Alpha')
axes[1].set_xlim(0, 90)
axes[1].legend(frameon=False)

plt.tight_layout()
plt.show()

# --- Norm of Cross Products ---
plt.figure(figsize=(5, 4))
plt.plot(alpha*180/np.pi, norm_cdot[:,0], label='|Cross 1-3|')
plt.plot(alpha*180/np.pi, norm_cdot[:,1], label='|Cross 1-4|')
plt.plot(alpha*180/np.pi, norm_cdot[:,2], label='|Cross 1-2|')
plt.xlabel(r'$\alpha$ (degrees)')
plt.ylabel('Norm of Cross Product')
plt.title('Norm of Cross Products vs Alpha')
plt.xlim(0, 90)
plt.legend(frameon=False)
plt.tight_layout()
plt.show()