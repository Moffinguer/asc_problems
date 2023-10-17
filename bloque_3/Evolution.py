#!/usr/bin/python3

import sys
import os
import random
import math
import matplotlib.pyplot as plt

def constraint_sigma_evolution(gen, growth_factor):
    return 1 / (1 + math.exp(-growth_factor * gen))

def constraints(problem, sol, population, generation):
    if problem == "zdt3":
        return 0

    restriction1 = sol[1] - 0.8 * sol[0] * math.sin(6 * math.pi * sol[0] + 2 * math.pi / len(population))
    restriction2 = sol[3] - 0.8 * sol[0] * math.sin(6 * math.pi * sol[0] + 4 * math.pi / len(population))

    sg1 = 0.5 * (1 - sol[0]) - (1 - sol[0]) ** 2
    sg2 = 0.25 * math.sqrt(1 - sol[0]) - 0.5 * (1 - sol[0])

    if sg1 != 0:
        restriction1 -= (sg1 / abs(sg1)) * math.sqrt(abs(0.5 * (1 - sol[0]) - (1 - sol[0]) ** 2))

    if sg2 != 0:
        restriction2 -= (sg2 / abs(sg2)) * math.sqrt(abs(0.25 * math.sqrt(1 - sol[0]) - 0.5 * (1 - sol[0])))

    restricted = 0
    if restriction1 < 0:
        restricted += (generation + 1) * constraint_sigma_evolution(generation, 0.1)

    if restriction2 < 0:
        restricted += (generation + 1) * constraint_sigma_evolution(generation, 0.2)

    return restricted

def mutate(x, window, upper_limit, lower_limit, dimensions, population, algorithm):
    CR, F, PR, SIG = 0.5, 0.5, 1 / len(population), 20
    chromosome_mutated = []
    parents = [window[random.randint(0, len(window) - 1)][2] for _ in range(3)]
    sigma = (upper_limit - lower_limit) / SIG

    for i in range(dimensions):
        # Cross
        if random.random() <= CR or i == int(random.random() * (dimensions - 1)):
            info = population[parents[0]][i] + F * (population[parents[1]][i] - population[parents[2]][i])
        else:
            info = x[i]

        # Gauss
        #info += math.exp(-0.5 * (info / sigma)**2) / (sigma * math.sqrt(2 * math.pi)) if random.random() <= PR

        chromosome_mutated.append(info)

    if algorithm == "zdt3":
        for i in range(len(chromosome_mutated)):
            while chromosome_mutated[i] < lower_limit or chromosome_mutated[i] > upper_limit:
                if chromosome_mutated[i] <= lower_limit:
                    chromosome_mutated[i] = lower_limit + (lower_limit - chromosome_mutated[i])
                if chromosome_mutated[i] >= upper_limit:
                    chromosome_mutated[i] = upper_limit - (chromosome_mutated[i] - upper_limit)
    else:
        while chromosome_mutated[0] < lower_limit or chromosome_mutated[0] > upper_limit:
            if chromosome_mutated[0] <= lower_limit:
                chromosome_mutated[0] = lower_limit + (lower_limit - chromosome_mutated[0])
            if chromosome_mutated[0] >= upper_limit:
                chromosome_mutated[0] = upper_limit - (chromosome_mutated[0] - upper_limit)

        cf6_lower_limit = -2
        cf6_upper_limit = 2
        for i in range(1, len(chromosome_mutated)):
            while chromosome_mutated[i] < cf6_lower_limit or chromosome_mutated[i] > cf6_upper_limit:
                if chromosome_mutated[i] <= cf6_lower_limit:
                    chromosome_mutated[i] = cf6_lower_limit + (cf6_lower_limit - chromosome_mutated[i])
                if chromosome_mutated[i] >= cf6_upper_limit:
                    chromosome_mutated[i] = cf6_upper_limit - (chromosome_mutated[i] - cf6_upper_limit)

    return chromosome_mutated

def select_subproblems(index, vecinity, population, lambda_window):
    number_near = int(population * vecinity)

    # print(f"SELECT {number_near} NEIGHBOURS NEAR {index} LAMBDA[{lambda_window[index][0]}, {lambda_window[index][1]}]\n")

    new_window = sorted(lambda_window, key=lambda a: euclidean_distance(lambda_window[index][0], a[0], lambda_window[index][1], a[1]))

    # print(new_window)

    new_window = new_window[:number_near]

    return new_window

def initialization(n, upper_limit, lower_limit, dimensions, algorithm):
    population_list = []
    chromosome = []

    if algorithm == "zdt3":
        for _ in range(n):
            for _ in range(dimensions):
                alelo = random.uniform(lower_limit, upper_limit)
                chromosome.append(alelo)
            population_list.append(chromosome)
            chromosome = []
    else:
        for _ in range(n):
            alelo = random.uniform(lower_limit, upper_limit)
            chromosome.append(alelo)
            for _ in range(1, dimensions):
                alelo = random.uniform(-2, 2)
                chromosome.append(alelo)
            population_list.append(chromosome)
            chromosome = []

    #print(f"GENERATED {n} INDIVIDUALS WITH {dimensions} CHROMOSOMES")

    return population_list

# Function zdt3
def zdt3(x, n):
    f_1 = x[0]
    g = sum(x[1:])
    g = 9 * g / (n - 1) + 1
    h = f_1 / g
    h = 1 - (h ** 0.5) - h * math.sin(10 * f_1 * math.pi)
    f_2 = g * h
    return f_1, f_2

# Function cf6
def cf6(x, n):
    f_1, f_2 = x[0], (1 - x[0]) ** 2
    for i in range(1, n - 1):
        if i % 2 == 1:
            f_1 += (x[i] - 0.8 * x[0] * math.cos(6 * math.pi * x[0] + i * math.pi / n)) ** 2
        else:
            f_2 += (x[i] - 0.8 * x[0] * math.sin(6 * math.pi * x[0] + i * math.pi / n)) ** 2
    return f_1, f_2

# Problem function
def problem_function(x, n, problem):
    if problem == "zdt3":
        return zdt3(x, n)
    else:
        return cf6(x, n)


def window(p):
    step = 1 / (p - 1)
    lambda_window = []
    lambd = 0
    for i in range(p):
        # Add a list of N values, where each subproblem is incremented by 1 with the index of the individual
        lambda_window.append([lambd, 1 - lambd, i])
        lambd += step
    return lambda_window

def euclidean_distance(x_00, x_01, x_10, x_11):
    dist = ((x_00 - x_01) ** 2 + (x_10 - x_11) ** 2) ** 0.5
    return dist

def read_input(filename):
    population, generations, neighborhood, inferior_limit, upper_limit, dimensions, experiments, algorithm = 100, 100, 0.25, 0, 1, 30, 1, "zdt3"

    filename = f"./INPUT_FILES/{filename}.in"
    if not os.path.exists(filename):
        print(f"{filename} does not exist")
        return population, generations, neighborhood, inferior_limit, upper_limit, dimensions, experiments, algorithm

    try:
        with open(filename, 'r') as input_file:
            for line in input_file:
                import re
                match = re.match(r'\s*(algorithm|experiments|population|generations|dimensions|upperLimit|inferiorLimit|neighborhood)=(\S+)\s*', line)
                if match:
                    key, value = match.groups()
                    if key == "population":
                        population = int(value)
                    elif key == "experiments":
                        experiments = int(value)
                    elif key == "generations":
                        generations = int(value)
                    elif key == "dimensions":
                        dimensions = int(value)
                    elif key == "upperLimit":
                        upper_limit = float(value)
                    elif key == "inferiorLimit":
                        inferior_limit = float(value)
                    elif key == "neighborhood":
                        neighborhood = float(value)
                    elif key == "algorithm":
                        algorithm = value
    except Exception as e:
        print(f"Error opening {filename}: {e}")
        sys.exit(1)

    return population, generations, neighborhood, inferior_limit, upper_limit, dimensions, experiments, algorithm


filename = sys.argv[1] if len(sys.argv) > 1 else "default_input"
population, generations, neighborhood, inferior_limit, upper_limit, dimensions, experiments, algorithm = read_input(filename)

print("####################################")
print("Variable\tValor")
print(f"ALGORITHM\t{algorithm}")
print(f"POPULATION\t{population}")
print(f"GENERATIONS\t{generations}")
print(f"NEIGHBORHOOD\t{neighborhood}")
print(f"INFERIOR_LIMIT\t{inferior_limit}")
print(f"UPPER_LIMIT\t{upper_limit}")
print(f"DIMENSIONS\t{dimensions}")
print(f"EXPERIMENTS\t{experiments}")
print("###################################")

lambda_window, population_list, f_1, f_2, evaluated_functions, z_1, z_2, times = (
        {}, [], None, None, None, None, None, None
)

directory_path = "./results/ZDT3"
if algorithm == "cf6":
    directory_path = "./results/CF6"

output_file = None


for experiment in range(experiments):

    plt.figure()
    plt.xlabel('f_1')
    plt.ylabel('f_2')
    plt.title('Gráfico de Puntos Dinámico')

    print("...................\nEXECUTION ", experiment, "...\n")

    if not os.path.exists(directory_path):
        try:
            os.makedirs(directory_path)
        except OSError as e:
            print("Unable to create", directory_path, ":", e)
            sys.exit(1)
    output_file = "./results/ZDT3/" + str(sys.argv[1]) + str(experiment) + ".out"
    if not os.path.exists(output_file):
        try:
            with open(output_file, 'w') as create_file:
                pass
        except IOError as e:
            print("Error opening", output_file, ":", e)
            sys.exit(1)
    else:
        try:
            with open(output_file, 'w') as clear_file:
                pass
        except IOError as e:
            print("Error cleaning", output_file, ":", e)
            sys.exit(1)

    lambda_value = window(population)
    population_list = initialization(population, upper_limit, inferior_limit, dimensions, algorithm)
    z_1, z_2 = None, None

    for generation in range(generations):

        print("...................\nGENERATION ", generation, "...\n")
        for x in population_list:
            f_1, f_2 = problem_function(x, dimensions, algorithm)
            violations = constraints(algorithm, x, population_list, generation)
            f_1 -= violations
            f_2 -= violations
            if z_1 is None or f_1 < z_1:
                z_1 = f_1
            if z_2 is None or f_2 < z_2:
                z_2 = f_2

		# Time to evolve
        for pop in range(population):
            if pop not in lambda_window:
                lambda_window[pop] = select_subproblems(pop, neighborhood, population, lambda_value)

            # Reproduce
            new_sols = mutate(population_list[pop], lambda_window[pop], upper_limit, inferior_limit, dimensions, population_list, algorithm)

            # Evaluation
            f_1, f_2 = problem_function(new_sols, dimensions, algorithm)
            violations = constraints(algorithm, new_sols, population_list, generation)
            f_1 -= violations
            f_2 -= violations

            with open(output_file, "a") as file_handle:
                file_handle.write(f"{f_1}\t{f_2}\t{violations}\n")

            # Update Best Sol
            z_1 = min(z_1, f_1)
            z_2 = min(z_2, f_2)

            # Update Neighbors
            # 0 -> [(A_0, A_1), (j[0], j[1]), (0, 1, 0), (0.5, 0.5, 1), (1, 0, 2)]
            for j in lambda_window[pop]:
                tchebycheff1 = j[0] * abs(f_1 - z_1)
                tch11 = j[1] * abs(f_2 - z_2)
                if tchebycheff1 < tch11:
                    tchebycheff1 = tch11

                f_1t, f_2t = problem_function(population_list[j[2]], dimensions, algorithm)

                tchebycheff2 = j[0] * abs(f_1t - z_1)
                tch22 = j[1] * abs(f_2t - z_2)
                if tchebycheff2 < tch22:
                    tchebycheff2 = tch22

                if tchebycheff1 > tchebycheff2:
                    population_list[j[2]] = new_sols

            tchebycheff1 = None
            if generation == generations - 1:
                plt.scatter(f_1, f_2, color='blue', marker='o', label='Puntos')

    lambda_window = {}
    population_list = []
    z_1, z_2 = None, None
    new_sols = None
    plt.show()

    f_1t, f_2t = None, None
