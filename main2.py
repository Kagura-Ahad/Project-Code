import itertools
import math
import time
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import csr_matrix
from networkx import Graph, is_connected, eulerian_circuit
import matplotlib.pyplot as plt
import numpy as np

class City:
    def __init__(self, name, coords):
        self.name = name
        self.coords = coords

# Function to calculate distance between two points
def calculate_distance(p1, p2):
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

# Function to calculate total distance of a path
def total_distance(path):
    return sum(calculate_distance(path[i - 1].coords, path[i].coords) for i in range(len(path)))

# Brute force TSP
def tsp_brute_force(cities):

    #Initisalise the minimum path as None so that we can check if it has been set
    min_path = None

    #Initialise the minimum distance as infinity so that we can check if it has been set
    min_distance = float('inf')

    #Start the timer
    start_time = time.time()

    #Iterate through all permutations of the cities
    for path in itertools.permutations(cities):
        current_distance = total_distance(path)
        if current_distance < min_distance:
            min_distance = current_distance
            min_path = path
    
    #Convert the path to a list
    min_path = list(min_path)

    #return to start
    min_path.append(min_path[0])

    #End the timer
    end_time = time.time()

    #Calculate the time taken
    time_taken = end_time - start_time

    #Return the minimum path, minimum distance and time taken
    return min_path, min_distance, time_taken

# Greedy TSP
def tsp_greedy(cities):
    start_city = cities[0]
    path = [start_city]
    unvisited_cities = set(cities[1:])
    current_city = start_city

    start_time = time.time()

    while unvisited_cities:
        next_city = min(unvisited_cities, key=lambda city: calculate_distance(current_city.coords, city.coords))
        unvisited_cities.remove(next_city)
        path.append(next_city)
        current_city = next_city

    path.append(start_city)  # return to start

    total_distance = sum(calculate_distance(city1.coords, city2.coords) for city1, city2 in zip(path[:-1], path[1:]))

    end_time = time.time()
    time_taken = end_time - start_time

    return path, total_distance, time_taken

# Dynamic programming TSP
def tsp_dynamic_programming(cities):
    n = len(cities)
    dist = [[calculate_distance(cities[i].coords, cities[j].coords) for j in range(n)] for i in range(n)]
    dp = [[None] * n for _ in range(1 << n)]
    path = [[0] * n for _ in range(1 << n)]

    def tsp(mask, pos):
        if mask == (1 << n) - 1:
            return dist[pos][0]
        if dp[mask][pos] is not None:
            return dp[mask][pos]
        ans = float('inf')
        for nxt in range(n):
            if mask >> nxt & 1:
                continue
            cur = dist[pos][nxt] + tsp(mask | (1 << nxt), nxt)
            if cur < ans:
                ans = cur
                path[mask][pos] = nxt
        dp[mask][pos] = ans
        return ans

    start_time = time.time()
    min_distance = tsp(1, 0)
    end_time = time.time()
    time_taken = end_time - start_time

    # Reconstruct the path
    mask = 1
    pos = 0
    min_path = [cities[0]]
    while len(min_path) < n:
        pos = path[mask][pos]
        mask |= 1 << pos
        min_path.append(cities[pos])
    min_path.append(cities[0])  # return to start

    return min_path, min_distance, time_taken

# Function to find a minimum weight matching for the odd degree vertices
def minimum_weight_matching(MST, G, odd_vert):
    while odd_vert:
        v = odd_vert.pop()
        distance = float('inf')
        u = 1
        closest = 0
        for u in odd_vert:
            if G[v][u] < distance:
                distance = G[v][u]
                closest = u
        MST.add_edge(v, closest, weight=distance)
        odd_vert.remove(closest)

# Function to find odd degree vertices in the graph
def find_odd_vertex(MST):
    odd = []
    for i in MST.nodes():
        if MST.degree(i) % 2 != 0:
            odd.append(i)
    return odd

# Christofides algorithm
def christofides(G, pos, cities):
    start_time = time.time()

    # Create a minimum spanning tree of graph G
    T = Graph(minimum_spanning_tree(csr_matrix(create_adjacency_matrix(G))))
    
    # Find all vertices of odd degree in the MST
    O = find_odd_vertex(T)
    
    # Add minimum weight matching edges to T
    minimum_weight_matching(T, G, O)
    
    # Create an Eulerian circuit
    euler_circuit = list(eulerian_circuit(T, source=pos))

    # Make the circuit found into a path
    path = [u for u, v in euler_circuit]
    path.append(path[0])

    end_time = time.time()
    time_taken = end_time - start_time

    # Calculate the total distance of the path
    total_distance = sum(G[path[i]][path[i + 1]] for i in range(len(path) - 1))

    #Convert the path to city objects
    path = [cities[city] for city in path]

    return path, total_distance, time_taken

# Function to create a graph from city coordinates
def create_graph(cities):
    G = {}
    for i, city1 in enumerate(cities):
        for j, city2 in enumerate(cities):
            if i != j:
                if i not in G:
                    G[i] = {}
                G[i][j] = calculate_distance(city1.coords, city2.coords)
    return G

# Assuming G is a dictionary of dictionaries representing the adjacency matrix
def create_adjacency_matrix(G):
    n = len(G)
    adjacency_matrix = [[0 if i == j else G[i].get(j, float('inf')) for j in range(n)] for i in range(n)]
    return adjacency_matrix

# Function to read data from a file
def read_data(file_name, max_cities):
    with open(file_name, 'r') as file:
        data = file.readlines()
    return [line.strip() for line in data[:max_cities]]

# Function to read coordinates from a file
def read_coords(file_name, max_cities):
    with open(file_name, 'r') as file:
        data = file.readlines()
    return [tuple(map(float, line.split())) for line in data[:max_cities]]

# Maximum number of cities to read
max_cities = 16  # Adjust this to 8, 9, or 10 for your desired number of cities

# Load your data from the files
city_names = read_data('names.txt', max_cities)
city_coords = read_coords('cyl_XY_coords.txt', max_cities)

# Create City objects
cities = [City(name, coords) for name, coords in zip(city_names, city_coords)]

# Call christofides algorithm for number of cities = 8
# G = create_graph(cities)
# path, totalDistance, time_taken = christofides(G, 0, cities)
# print(f"Minimum path: {[city.name for city in path]}\nMinimum distance: {totalDistance}\nTime taken: {time_taken} seconds")

# Function to run all algorithms and plot their time complexities
def run_algorithms_and_plot(cities):
    # Prepare data structures for plotting
    city_counts = range(4, 10)
    # brute_force_times = []
    # brute_force_distances = []
    # greedy_times = []
    # greedy_distances = []
    # dynamic_times = []
    # dynamic_distances = []
    christofides_times = []
    christofides_distances = []

    # Run each algorithm for different city counts and record the time taken
    for count in city_counts:
        selected_cities = cities[:count]

        # _, brute_force_distance, brute_force_time = tsp_brute_force(selected_cities)
        # brute_force_distances.append(brute_force_distance)
        # brute_force_times.append(brute_force_time)

        # _, greedy_distance, greedy_time = tsp_greedy(selected_cities)
        # greedy_times.append(greedy_time)
        # greedy_distances.append(greedy_distance)

        # _, dynamic_distance, dynamic_time = tsp_dynamic_programming(selected_cities)
        # dynamic_times.append(dynamic_time)
        # dynamic_distances.append(dynamic_distance)

        G = create_graph(selected_cities)
        path, total_distance, time_taken = christofides(G, 0, selected_cities)
        christofides_times.append(time_taken)
        christofides_distances.append(total_distance)

    plt.figure(figsize=(14, 8))
    # Plotting Time Complexity
    plt.subplot(1, 2, 1)
    # plt.plot(city_counts, brute_force_times, label='Brute Force', marker='o')
    # plt.plot(city_counts, greedy_times, label='Greedy', marker='s')
    # plt.plot(city_counts, dynamic_times, label='Dynamic Programming', marker='^')
    plt.plot(city_counts, christofides_times, label='Christofides', marker='x')
    plt.xlabel('Number of Cities')
    plt.ylabel('Time Taken (seconds)')
    plt.title('TSP Algorithms Time Complexity')
    plt.legend()
    plt.grid(True)

    # Plotting Total Distance
    plt.subplot(1, 2, 2)
    # plt.plot(city_counts, brute_force_distances, label='Brute Force', marker='o')
    # plt.plot(city_counts, greedy_distances, label='Greedy', marker='s')
    # plt.plot(city_counts, dynamic_distances, label='Dynamic Programming', marker='^')
    plt.plot(city_counts, christofides_distances, label='Christofides', marker='x')
    plt.xlabel('Number of Cities')
    plt.ylabel('Total Distance')
    plt.title('TSP Algorithms Total Distance')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()  # Adjust spacing between subplots
    plt.show()

# Run the function
run_algorithms_and_plot(cities)