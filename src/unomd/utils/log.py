import matplotlib.pyplot as plt
import numpy as np

# Function to parse the log file
def parse_log_file(log_file_path):
    print("Parsing log file")
    steps = []
    energies = []
    temperatures = []  # List to store temperature values
    with open(log_file_path, 'r') as file:
        # Skip header or lines without numeric data
        next(file)
        for line in file:
            cols = line.strip().split('\t')
            # print(cols[2])
            try:
                steps.append(float(cols[2]))
                energies.append(float(cols[3]))
                temperatures.append(float(cols[4]))  # Parse temperature
            except ValueError as e:
                print(e)
    return steps, energies, temperatures

def plot_energy_temp(steps, energies, temperatures):
    # Create figure and first axis
    plt.figure(figsize=(10, 6))
    ax1 = plt.gca()  # Get current axis
    ax2 = ax1.twinx()  # Create another axis that shares the same x-axis
    
    # Plot energy on the first y-axis
    ax1.plot(steps, energies, marker='o', linestyle='-', color='blue', label='Potential Energy')
    ax1.set_xlabel('Time Step')
    ax1.set_ylabel('Energy (kJ/mol)', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    
    # Plot temperature on the second y-axis
    ax2.plot(steps, temperatures, marker='x', linestyle='-', color='red', label='Temperature')
    ax2.set_ylabel('Temperature (F)', color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    # Flip the temperature axis
    ax1.invert_yaxis()
    
    # Title and grid
    plt.title('Energy and Temperature vs. Time Step')
    ax1.grid(True)
    
    # Optional: add a legend. Comment these lines if you find the legend unnecessary.
    ax1.legend(loc='upper left')
    # ax2.legend(loc='upper right')
    
    plt.show()

def analyze_distributions(temperatures, energies):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # Temperature distribution
    ax1.hist(temperatures, bins=50, density=True, alpha=0.7)
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Density')
    ax1.set_title('Temperature Distribution')
    
    # Energy distribution
    ax2.hist(energies, bins=50, density=True, alpha=0.7)
    ax2.set_xlabel('Potential Energy (kJ/mol)')
    ax2.set_ylabel('Density')
    ax2.set_title('Energy Distribution')
    
    plt.tight_layout()
    plt.show()

def calculate_running_averages(times, temperatures, energies, window=10):
    temp_running_avg = []
    energy_running_avg = []
    
    for i in range(len(times)):
        start_idx = max(0, i - window)
        temp_running_avg.append(np.mean(temperatures[start_idx:i+1]))
        energy_running_avg.append(np.mean(energies[start_idx:i+1]))
    
    plt.figure(figsize=(10, 6))
    plt.plot(times, temp_running_avg, label='Temperature Running Average')
    plt.xlabel('Time (ps)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.grid(True)
    plt.show()

def calculate_autocorrelation(data):
    data = np.array(data)
    data = data - np.mean(data)
    autocorr = np.correlate(data, data, mode='full')
    autocorr = autocorr[len(autocorr)//2:]
    autocorr = autocorr / autocorr[0]
    
    # Calculate decorrelation time
    decay_idx = np.where(autocorr < np.exp(-1))[0]
    if len(decay_idx) > 0:
        decorr_time = decay_idx[0]
    else:
        decorr_time = None
    
    return autocorr, decorr_time

def plot_autocorrelation(temperatures, energies):
    temp_autocorr, temp_decorr = calculate_autocorrelation(temperatures)
    energy_autocorr, energy_decorr = calculate_autocorrelation(energies)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    ax1.plot(temp_autocorr[:1000])
    ax1.set_title('Temperature Autocorrelation')
    if temp_decorr:
        ax1.axvline(x=temp_decorr, color='r', linestyle='--', 
                    label=f'Decorrelation time: {temp_decorr} steps')
    ax1.legend()
    
    ax2.plot(energy_autocorr[:1000])
    ax2.set_title('Energy Autocorrelation')
    if energy_decorr:
        ax2.axvline(x=energy_decorr, color='r', linestyle='--',
                    label=f'Decorrelation time: {energy_decorr} steps')
    ax2.legend()
    
    plt.tight_layout()
    plt.show()