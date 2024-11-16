import numpy as np
import re
import glob

def read_zone_file(filename):
    data = []
    with open(filename, 'r') as file:
        content = file.readlines()

    iteration_data = {}
    i = 0
    while i < len(content):
        line = content[i].strip()
        if line.startswith('Iteration:'):
            # Start a new iteration data block
            if iteration_data:
                data.append(iteration_data)
            iteration_data = {}
            iteration_number = int(line.split(':')[1].strip())
            iteration_data['Iteration'] = iteration_number
            i += 1
        elif ':' in line:
            key = line.split(':')[0].strip()
            if key in ['eps_nuc', 'rho', 'T', 'd_eps_nuc_dRho', 'd_eps_nuc_dT']:
                # Scalar variable
                value = float(line.split(':')[1].strip())
                iteration_data[key] = value
                i += 1
            elif key in ['d_epsnuc_dx', 'dxdt_nuc', 'd_dxdt_nuc_dRho', 'd_dxdt_nuc_dT']:
                # Vector variable
                i += 1  # Move to the next line where data starts
                values = []
                while i < len(content) and not content[i].strip().endswith(':'):
                    line_values = content[i].strip().split()
                    values.extend([float(val) for val in line_values])
                    i += 1
                iteration_data[key] = np.array(values)
            elif key == 'd_dxdt_nuc_dx':
                # 2D array variable
                i += 1  # Move to the next line where data starts
                rows = []
                while i < len(content) and not content[i].strip().startswith(('Iteration:', 'eps_nuc:', 'rho:', 'T:', 'd_eps_nuc_dRho:', 'd_eps_nuc_dT:', 'd_epsnuc_dx:', 'dxdt_nuc:', 'd_dxdt_nuc_dRho:', 
'd_dxdt_nuc_dT:')):
                    line_values = content[i].strip().split()
                    if line_values:
                        row = [float(val) for val in line_values]
                        rows.append(row)
                    i += 1
                iteration_data[key] = np.array(rows)
            else:
                # Unknown key
                i += 1
        else:
            i += 1

    # Append the last iteration data
    if iteration_data:
        data.append(iteration_data)

    return data

def read_all_zones():
    zone_files = sorted(glob.glob('zone_*.txt'))
    all_zone_data = {}
    for filename in zone_files:
        zone_number = re.findall(r'zone_(\d+)\.txt', filename)[0]
        zone_data = read_zone_file(filename)
        all_zone_data[zone_number] = zone_data
    return all_zone_data

# Example usage
if __name__ == '__main__':
    all_data = read_all_zones()

    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 6))
    for i in [1815,1816,1817,1818,1819]:
        # Access data for a specific zone and iteration
        zone_number = str(i)  # Change as needed
        zone_data = all_data[zone_number]

        # For example, get the eps_nuc values over iterations
        iterations = [iter_data['Iteration'] for iter_data in zone_data]
        eps_nuc_values = [iter_data['eps_nuc'] for iter_data in zone_data]

        # Plotting eps_nuc over iterations for the specified zone

        eps_nuc_values = np.array(eps_nuc_values)
        # Initialize an array to store the log-scaled values
        eps_nuc_values_log = np.zeros_like(eps_nuc_values)

        # Handle positive terms: directly take the logarithm
        positive_mask = eps_nuc_values > 0
        eps_nuc_values_log[positive_mask] = np.log10(eps_nuc_values[positive_mask])

        # Handle negative terms: take the logarithm of the absolute value and scale by -1
        negative_mask = eps_nuc_values < 0
        eps_nuc_values_log[negative_mask] = -np.log10(np.abs(eps_nuc_values[negative_mask]))
        
        #name = "zone" +str(i)
        plt.plot(iterations,eps_nuc_values_log, marker='o', label = str(i))
    

#    plt.plot(iterations,eps_nuc_values_log, marker='o')
#    plt.ylim(-1e21,1e21)
    print (iterations)
#    print (eps_nuc_values_log)
    plt.xlabel('Iteration')
    plt.ylabel('eps_nuc')
    plt.title(f'eps_nuc over all iterations for core zones')
    plt.legend()
    plt.grid(True)
    plt.savefig('plots/eps_nuc.pdf',dpi = 300)
    plt.show()
    





    "plot partials for d_eps_nuc_dT in core zones"
    plt.figure(figsize=(8, 6))
    for i in [1815,1816,1817,1818,1819]:
        # Access data for a specific zone and iteration
        zone_number = str(i)  # Change as needed
        zone_data = all_data[zone_number]

        # For example, get the eps_nuc values over iterations
        iterations = [iter_data['Iteration'] for iter_data in zone_data]
        eps_nuc_values = [iter_data['d_eps_nuc_dT'] for iter_data in zone_data]

        # Plotting eps_nuc over iterations for the specified zone

        eps_nuc_values = np.array(eps_nuc_values)
        # Initialize an array to store the log-scaled values
        eps_nuc_values_log = np.zeros_like(eps_nuc_values)

        # Handle positive terms: directly take the logarithm
        positive_mask = eps_nuc_values > 0
        eps_nuc_values_log[positive_mask] = np.log10(eps_nuc_values[positive_mask])

        # Handle negative terms: take the logarithm of the absolute value and scale by -1
        negative_mask = eps_nuc_values < 0
        eps_nuc_values_log[negative_mask] = -np.log10(np.abs(eps_nuc_values[negative_mask]))
        
        #name = "zone" +str(i)
        plt.plot(iterations,eps_nuc_values_log, marker='o', label = str(i))
    

#    plt.plot(iterations,eps_nuc_values_log, marker='o')
#    plt.ylim(-1e21,1e21)
    print (iterations)
#    print (eps_nuc_values_log)
    plt.xlabel('Iteration')
    plt.ylabel('d_eps_nuc_dT')
    plt.title(f'd_eps_nuc_dT over all iterations for core zones')
    plt.legend()
    plt.grid(True)
    plt.savefig('plots/d_eps_nuc_dT.pdf',dpi = 300)
    plt.show()
    


    "plot partials for d_eps_nuc_dRho in core zones"
    plt.figure(figsize=(8, 6))
    for i in [1815,1816,1817,1818,1819]:
        # Access data for a specific zone and iteration
        zone_number = str(i)  # Change as needed
        zone_data = all_data[zone_number]

        # For example, get the eps_nuc values over iterations
        iterations = [iter_data['Iteration'] for iter_data in zone_data]
        eps_nuc_values = [iter_data['d_eps_nuc_dRho'] for iter_data in zone_data]

        # Plotting eps_nuc over iterations for the specified zone

        eps_nuc_values = np.array(eps_nuc_values)
        # Initialize an array to store the log-scaled values
        eps_nuc_values_log = np.zeros_like(eps_nuc_values)

        # Handle positive terms: directly take the logarithm
        positive_mask = eps_nuc_values > 0
        eps_nuc_values_log[positive_mask] = np.log10(eps_nuc_values[positive_mask])

        # Handle negative terms: take the logarithm of the absolute value and scale by -1
        negative_mask = eps_nuc_values < 0
        eps_nuc_values_log[negative_mask] = -np.log10(np.abs(eps_nuc_values[negative_mask]))
        
        #name = "zone" +str(i)
        plt.plot(iterations,eps_nuc_values_log, marker='o', label = str(i))
    

#    plt.plot(iterations,eps_nuc_values_log, marker='o')
#    plt.ylim(-1e21,1e21)
    print (iterations)
#    print (eps_nuc_values_log)
    plt.xlabel('Iteration')
    plt.ylabel('d_eps_nuc_dRho')
    plt.title(f'd_eps_nuc_dRho over all iterations for core zones')
    plt.legend()
    plt.grid(True)
    plt.savefig('plots/d_eps_nuc_dRho.pdf',dpi = 300)
    plt.show()
    
