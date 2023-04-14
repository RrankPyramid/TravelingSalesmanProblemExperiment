import subprocess
import os
import matplotlib.pyplot as plt
import sys

def compile_cpp_files(file_names):
    for file_name in file_names:
        result = subprocess.run(["g++", "-O3", "-o", file_name[:-4], file_name], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Failed to compile {file_name}. Error message:")
            print(result.stderr)
            sys.exit(1)

def read_dataset(file_path):
    cities = []
    with open(file_path, 'r') as file:
        for line in file.readlines():
            if line.strip() == 'EOF':
                continue
            index, x, y = map(float, line.strip().split())
            cities.append((x, y))
    return cities

def write_to_file(file_path, string):
    with open(file_path, 'w') as file:
        file.write(string)

def run_cpp_implementation(implementation_path, cities):
    input_data = f"{len(cities)}\n"
    for city in cities:
        input_data += f"{city[0]} {city[1]}\n"
#    write_to_file("./testdata/Djibouti/dj38.in", input_data)
#   exit(0)
    result = subprocess.run([implementation_path], input=input_data, text=True, capture_output=True)
    return result.stdout

def test(dataset_path, dataset_name):
    implementations = [("Tabu", "./Tabu"), ("ACO", "./ACO"), ("Hybrid", "./Hybrid")]
    cities = read_dataset(dataset_path)
    
    result = {}
    for name, implementation in implementations:
        output = run_cpp_implementation(implementation, cities)
        print("Output of", name, "implementation:", output)
        result_length = float(output.split("Best tour length: ")[1].split("\n")[0])
        print("Best tour length:", result_length)
        result[name] = result_length

    return dataset_name, result

import numpy as np
def show(running_result_list):
    fig, ax = plt.subplots()
    dataset_names = [entry[0] for entry in running_result_list]
    num_datasets = len(dataset_names)
    num_algorithms = len(running_result_list[0][1])
    bar_width = 0.2

    # Compute bar positions
    ind = np.arange(num_datasets)
    bar_positions = [ind + bar_width * i for i in range(num_algorithms)]

    # Plot bars for each algorithm
    for i, name in enumerate(running_result_list[0][1].keys()):
        results = [entry[1][name] for entry in running_result_list]
        print("Results for", name, ":", results)
        ax.bar(bar_positions[i], results, width=bar_width, label=name)

    ax.set_xticks(ind + bar_width * (num_algorithms - 1) / 2)
    ax.set_xticklabels(dataset_names)
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Best Tour Length")
    ax.set_title("Performance Comparison")
    ax.legend()

    plt.show()



if __name__ == "__main__":
    cpp_files = ["Tabu.cpp", "ACO.cpp", "Hybrid.cpp"]
    compile_cpp_files(cpp_files)

    datasets = [("./testdata/Djibouti/dj38.tsp", "Djibouti"), 
            ("./testdata/Luxembourg/lu980.tsp", "Luxembourg"), 
            ("./testdata/Qatar/qa194.tsp", "Qatar"), 
            ("./testdata/WesternSahara/wi29.tsp", "Western Sahara"),
            ("./testdata/Zimbabwe/zi929.tsp", "Zimbabwe"),
            ("./testdata/Uruguay/uy734.tsp", "Uruguay")
            ]    
    running_result_list = [test(dataset_path, dataset_name) for dataset_path, dataset_name in datasets]
    show(running_result_list)
