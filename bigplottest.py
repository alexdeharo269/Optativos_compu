import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Define the column names
columns = [
    "h", "vrmod", "radmult",
    "numPlanetas", "interacciones", "absorbidos", "supervivientes", "n_old", "numPlanetas_per_supervivientes",
    "smalls", "old_small", "avg_small_radius", "avg_small_periodo", "avg_small_exc", "avg_small_avgdtS", "avg_small_U", "avg_small_rock",
    "mediums", "old_mediums", "avg_medium_radius", "avg_medium_periodo", "avg_medium_exc", "avg_medium_avgdtS", "avg_medium_U", "avg_medium_rock",
    "larges", "old_larges", "avg_large_radius", "avg_large_periodo", "avg_large_exc", "avg_large_avgdtS", "avg_large_U", "avg_large_rock",
    "very_larges", "old_very_larges", "avg_very_large_radius", "avg_very_large_periodo", "avg_very_large_exc", "avg_very_large_avgdtS", "avg_very_large_U", "avg_very_large_rock"
]

# Read the multiverse file into a DataFrame
file_path = r'C:\Users\Ale\Desktop\resumen.dat'  # Adjust this path as needed
df = pd.read_csv(file_path, sep='\s+', header=None, names=columns)

# Function to plot multiple vrmod values on the same plot
def plot_multiple_vrmod(df, x_col, y_cols, filter_col, vrmod_values):
    plt.figure(figsize=(4, 2.5))
    ax = plt.gca()
    for vrmod in vrmod_values:
        filtered_df = df[df[filter_col] == vrmod]
        print(filtered_df['n_old'])
        for y_col in y_cols:
            plt.plot(filtered_df[x_col], filtered_df[y_col], marker='o', label=f'$v_r$={vrmod}')
    plt.xlabel('Radio', fontsize=13, fontname='Times New Roman')
    plt.ylabel("Distancia media al Sol", fontsize=10, fontname='Times New Roman')
    plt.title('Planetas muy grandes estables', fontsize=12, fontname='Times New Roman')
    plt.legend(loc='best', prop={'family': 'Times New Roman', 'size': 9})
    plt.grid(True)
    plt.xscale('log')
    ax.set_xlim([900, 60000])

    # Define a custom formatter for the y-axis to remove precision
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x)}'))

    plt.show()

# Example of plotting 'h' against multiple y-axis columns for specific vrmod values
vrmod_values = [0.1, 0.5, 1.0,5.0, 10.0]  # Specify the vrmod values you want to plot
y_columns = ["avg_very_large_avgdtS"]
plot_multiple_vrmod(df, 'radmult', y_columns, 'vrmod', vrmod_values)
