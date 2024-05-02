import math

NUM_CITIES = 4

class City:
    def __init__(self, name, x, y):
        self.name = name
        self.x = x
        self.y = y

def distance(a, b):
    return math.sqrt((a.x - b.x)**2 + (a.y - b.y)**2)

def swap(cities, i, j):
    cities[i], cities[j] = cities[j], cities[i]

def write_log(filename, algorithm, path, path_length, total_distance):
    with open(filename, "a") as file:
        file.write(f"Algorithm: {algorithm}\n")
        file.write(f"Total distance: {total_distance}\n")
        file.write("Path: ")
        for city in path:
            file.write(f"{city.name} -> ")
        file.write("\n\n")

def bruteForce(cities, l, r, best_path, best_distance):
    if l == r:
        total_distance = 0
        for i in range(r):
            total_distance += distance(cities[i], cities[i + 1])
        total_distance += distance(cities[r], cities[0]) # return to start
        if total_distance < best_distance[0]:
            best_distance[0] = total_distance
            best_path.clear()
            best_path.extend(cities[:r + 1])
            write_log("log.txt", "Brute Force", best_path, r + 1, total_distance)
    else:
        for i in range(l, r + 1):
            swap(cities[l], cities[i])
            bruteForce(cities, l + 1, r, best_path, best_distance)
            swap(cities[l], cities[i]) # backtrack

def greedy(cities, best_path, best_distance):
    visited = [False] * NUM_CITIES
    best_path.append(cities[0])
    visited[0] = True

    for i in range(1, NUM_CITIES):
        min_distance = float('inf')
        min_index = -1

        for j in range(NUM_CITIES):
            if not visited[j]:
                dist = distance(best_path[i - 1], cities[j])
                if dist < min_distance:
                    min_distance = dist
                    min_index = j

        best_path.append(cities[min_index])
        visited[min_index] = True
        best_distance[0] += min_distance
        write_log("log.txt", "Greedy", best_path, i + 1, best_distance[0])

    best_distance[0] += distance(best_path[NUM_CITIES - 1], best_path[0])
    write_log("log.txt", "Greedy", best_path, NUM_CITIES, best_distance[0])

def tsp(cities, best_path, best_distance, algorithm):
    if algorithm == "bruteForce":
        bruteForce(cities, 0, NUM_CITIES - 1, best_path, best_distance)
    elif algorithm == "greedy":
        greedy(cities, best_path, best_distance)

if __name__ == "__main__":
    cities = []
    best_path = []
    best_distance = [float('inf')]

    with open("names.txt", "r") as names_file, open("cyl_XY_coords.txt", "r") as coords_file:
        for _ in range(NUM_CITIES):
            name, _ = names_file.readline().strip().split(',')
            x, y = map(float, coords_file.readline().strip().split())
            cities.append(City(name, x, y))

    filename = "log.txt"
    with open(filename, "w") as file:
        file.write(f"Dataset contains {NUM_CITIES} cities\n\n")

    import time
    start = time.time()

    tsp(cities, best_path, best_distance, "greedy")

    end = time.time()

    time_taken = end - start

    with open("log.txt", "a") as file:
        file.write(f"Time taken: {time_taken}\n\n")
        file.write("Unit of time taken: seconds\n\n")

    print("Best path:", [city.name for city in best_path])
    print("Total distance:", best_distance[0])
