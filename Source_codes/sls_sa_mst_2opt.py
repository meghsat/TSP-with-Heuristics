import numpy as np
import random
import math
import time
import os
import glob
import csv

def calculate_total_distance(tour, distance_matrix):
    total_distance = 0
    for i in range(len(tour) - 1):
        total_distance += distance_matrix[tour[i]][tour[i + 1]]
    total_distance += distance_matrix[tour[-1]][tour[0]]  # Return to the starting city
    return total_distance

def distance(city1, city2):
    # Calculate Euclidean distance
    x_diff = city1[0] - city2[0]
    y_diff = city1[1] - city2[1]
    return math.sqrt(x_diff**2 + y_diff**2)

def generate_mst_solution(cities,num_cities):
    from collections import defaultdict
    start_node = np.argmin(np.min(cities, axis=1))
    # Create graph
    graph = defaultdict(list)
    for i, city1 in enumerate(cities):
        for j, city2 in enumerate(cities):
            if i != j:
                graph[i].append((j, distance(city1, city2)))

    # Find MST using Prim's algorithm
    visited = set()
    mst = set()
    current_node = start_node
    #print(current_node)
    visited.add(current_node)

    while len(visited) < len(cities):
        for neighbor, weight in graph[current_node]:
            if neighbor not in visited:
                min_edge = (current_node, neighbor, weight)
                break
        for neighbor, weight in graph[current_node]:
            if neighbor not in visited and weight < min_edge[2]:
                min_edge = (current_node, neighbor, weight)
        visited.add(min_edge[1])
        mst.add(min_edge)
        current_node = min_edge[1]

    # Extract tour from MST
    tour = [0]
    visited = set([0])
    while len(visited) < len(cities):
        for edge in mst:
            if edge[0] in visited and edge[1] not in visited:
                tour.append(edge[1])
                visited.add(edge[1])
                break

    return tour

def generate_initial_solution(num_cities):
    return random.sample(range(num_cities), num_cities)

def generate_neighbor_solution(current_solution):
    i, j = sorted(random.sample(range(len(current_solution)), 2))
    new_solution = current_solution[:i] + list(reversed(current_solution[i:j + 1])) + current_solution[j + 1:]
    return new_solution

def acceptance_probability(old_distance, new_distance, temperature):
    if new_distance < old_distance:
        return 1.0
    return math.exp((old_distance - new_distance) / (temperature + 1e-10))

def simulated_annealing(distance_matrix,num_cities, max_iterations=10000, initial_temperature=100, cooling_rate=0.9):
    #num_cities = len(distance_matrix)
    #current_solution = generate_initial_solution(num_cities)
    current_solution = generate_mst_solution(distance_matrix,num_cities) 
    current_distance = calculate_total_distance(current_solution, distance_matrix)

    best_solution = current_solution.copy()
    best_distance = current_distance

    temperature = initial_temperature

    for iteration in range(max_iterations):
        neighbor_solution = generate_neighbor_solution(current_solution)
        neighbor_distance = calculate_total_distance(neighbor_solution, distance_matrix)

        probability = acceptance_probability(current_distance, neighbor_distance, temperature)

        if random.uniform(0, 1) < probability:
            current_solution = neighbor_solution
            current_distance = neighbor_distance

        if current_distance < best_distance:
            best_solution = current_solution.copy()
            best_distance = current_distance

        temperature *= cooling_rate

    return best_solution, best_distance


folder_path = "C:/Users/namra/Downloads/Competion/Competion"
files_list = os.listdir(folder_path)
city_files = glob.glob(os.path.join(folder_path, '*.txt'))
output_csv_file = "C:/Users/namra/Downloads/Competion/output.csv"
with open(output_csv_file, 'w', newline='') as csvfile:
    for filepath in city_files:
        #filepath = "C:/Users/namra/Downloads/Competion/Competion/tsp-problem-1000-200000-100-25-1.txt"
        with open(filepath, 'r') as file:
            num_cities = int(file.readline().strip())
            distances = np.loadtxt(file, skiprows=0)
        #print(distances)
        start_time = time.time()
        best_tour, best_distance = simulated_annealing(distances,num_cities)
        end_time = time.time()
        exe_time = end_time - start_time
        print("Best Tour:", best_tour)
        print("Best Distance:", best_distance)
        print("Time: ", exe_time)
    
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow([filepath, best_tour, best_distance, exe_time])