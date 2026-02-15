import matplotlib.pyplot as plt
import pandas as pd

trajectories = pd.read_csv("data/trajectory_fixed_20260212_214002.csv", 
                           header=0, 
                           index_col=None)

X_vals = list(trajectories['time'])
Y1_vals = list(trajectories['prey'])
Y2_vals = list(trajectories['predator'])

plt.figure(figsize=(10, 6))
plt.plot(X_vals, Y1_vals, color='blue', marker='o', label='Prey')
plt.plot(X_vals, Y2_vals, color='black', marker='s', label='Predator')

plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Prey and Predator Dynamics')
plt.legend()
plt.xlim(min(X_vals), max(X_vals))
plt.ylim(0, max(max(Y1_vals), max(Y2_vals)) * 1.1)
plt.tight_layout()
plt.show()