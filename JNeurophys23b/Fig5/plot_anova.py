import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv('C:/Users/cclab/Documents/MATLAB/layer_WM/tuning/tbl_for_anova.csv')

# Create the plot
plt.figure(figsize=(10, 6))
sns.pointplot(data=data, x='group', y='rate', hue='stage', ci=95, capsize=.1)

# Set the title and labels
plt.title('Interaction Between Group and Stage on Firing Rate')
plt.xlabel('Group')
plt.ylabel('Firing Rate')

# Show the plot
plt.show()
