import numpy as np
import torch
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from scipy.stats import binned_statistic_dd


print("reading data")
data = torch.load("/home/modelflows/DATASETS/JHC_HM3_LES_dataset_T_YCO_YCO2_YH2O_tensor.pt", weights_only=True).numpy()
data_coordinates = torch.load("/home/modelflows/DATASETS/cell_centres.pt",  weights_only=True).numpy() 
print("Done! \n")




temperature  = data[0, ...]
del data



x = data_coordinates[:, 0]
y = data_coordinates[:, 1]
z = data_coordinates[:, 2]

r = np.sqrt(x**2 + y**2)
theta = np.arctan2(y, x)
cylindrical_coords = np.stack([r, theta, z], axis=1)
print(cylindrical_coords.shape)

r_max = r.max()/5*2 # radial limit
#
z_min, z_max = 0,  data_coordinates[:,2].max()/5*4  # axial limits
theta_min, theta_max = -np.pi, np.pi  # full circle
mask = (r < r_max) & (z >= z_min) & (z <= z_max) & (theta >= theta_min) & (theta <= theta_max)

filtered_coords = cylindrical_coords[mask]
filtered_data = temperature[mask, :]



def check_for_nans(n_r, n_theta, n_z, coord_bins, filtered_data, r_max, theta_min, theta_max, z_min, z_max):
    r_edges = np.linspace(0, r_max, n_r)
    theta_edges = np.linspace(theta_min, theta_max, n_theta)
    z_edges = np.linspace(z_min, z_max, n_z)

    tensor_data = np.zeros((n_r - 1, n_theta - 1, n_z - 1, filtered_data.shape[1]))

    for t in range(filtered_data.shape[1]):
        binned, _, _ = binned_statistic_dd(
            coord_bins,
            filtered_data[:, t],
            statistic='mean',
            bins=[r_edges, theta_edges, z_edges]
        )
        if np.isnan(binned).any():
            return True  # Contains NaNs
        tensor_data[..., t] = binned

    return False  # Safe

# Initial values
n_r, n_theta, n_z = 40, 77, 160 #(already optimal)
step = 2
coord_bins = np.stack([r[mask], theta[mask], z[mask]], axis=1)
print("Starting optimization:")
print("Optimizing n_r...")
while True:
    if check_for_nans(n_r + step, n_theta, n_z, coord_bins, filtered_data, r_max, theta_min, theta_max, z_min, z_max):
        break
    n_r += step
best_nr = n_r
print(f"Best n_r: {n_r}")

print("Optimizing n_theta...")
while True:
    if check_for_nans(n_r, n_theta + step, n_z, coord_bins, filtered_data, r_max, theta_min, theta_max, z_min, z_max):
        break
    n_theta += step
best_ntheta = n_theta
print(f"Best n_theta: {n_theta}")

print("Optimizing n_z...")
while True:
    if check_for_nans(n_r, n_theta, n_z + step, coord_bins, filtered_data, r_max, theta_min, theta_max, z_min, z_max):
        break
    n_z += step
best_nz = n_z
print(f"Best n_z: {n_z}")

best_nr, best_ntheta, best_nz = n_r, n_theta, n_z

r_edges = np.linspace(0, r_max, best_nr)
theta_edges = np.linspace(theta_min, theta_max, best_ntheta)
z_edges = np.linspace(z_min, z_max, best_nz)

coord_bins = np.stack([r[mask], theta[mask], z[mask]], axis=1)

tensor_data = np.zeros((n_r - 1, n_theta - 1, n_z - 1, filtered_data.shape[1]))

for t in tqdm(range(filtered_data.shape[1])):
    binned, _, _ = binned_statistic_dd(
        coord_bins,
        filtered_data[:, t],
        statistic='mean', 
        bins=[r_edges, theta_edges, z_edges]
    )
    tensor_data[..., t] = binned

np.save("binned_temperature.npy", tensor_data)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].imshow(tensor_data[0, :, :, 0].T, cmap='inferno')
axes[0].set_title("Constant r")
axes[0].axis('off')

axes[1].imshow(tensor_data[:, 0, :, 0].T, cmap='inferno')
axes[1].set_title("Constant theta")
axes[1].axis('off')

axes[2].imshow(tensor_data[:, :, 0, 0].T, cmap='inferno')
axes[2].set_title("Constant z")
axes[2].axis('off')

plt.tight_layout()
plt.show()


