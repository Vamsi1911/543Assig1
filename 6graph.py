import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

def plot_and_save_csv(filename):
    # Read the CSV file
    df = pd.read_csv(filename)
    
    # Plot the data
    plt.figure(figsize=(10, 6))
    for column in df.columns[1:]:
        plt.plot(df['x'], df[column], label=column)
    
    plt.xlabel('Position (m)')
    plt.ylabel('Temperature (°C)')
    plt.title(f'Temperature Distribution - {os.path.basename(filename)}')
    plt.legend()
    plt.grid(True)
    
    # Save the plot as a PNG file
    output_filename = os.path.splitext(filename)[0] + '.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved as {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide a CSV filename as an argument.")
    else:
        plot_and_save_csv(sys.argv[1])