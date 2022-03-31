import copy
import math
import random
import matplotlib.pyplot as plt
import numpy as np

######input
dim_population = 0
p_mutation = 0
p_cross = 0
A = 0
B = 0
precision = 0
a = 0
b = 0
c = 0
steps = 0

######variables
population = []
lg = 0
intermediate_population = []    #after the selection
recombinated_population = []    #after the recombination
mutated_population = []         #after mutation
values_fittest = []     #fit of the fittest element
fittest_from_current_generation = []

is_elitist = 0 # 1 - the fittest elem is qualified for the next generation, 0 - is not
type_mutation = 0 # 0 - rare muation, 1 - normal mutation
ct = 0

def read():
    global dim_population, p_mutation, p_cross
    global A, B, precision #domenion of the function and the precision of the numbers
    global a, b, c #polynomial coefficients
    global steps #number of times that we have to run the algorithm

    f = open("input.txt", "r")
    dim_population = int(f.readline())
    A = float(f.readline())
    B = float(f.readline())
    a = float(f.readline())
    b = float(f.readline())
    c = float(f.readline())
    precision = int(f.readline())
    p_cross = float(f.readline())
    p_mutation = float(f.readline())
    steps = int(f.readline())
    f.close()

    global is_elitist, type_mutation
    is_elitist = int(input("Elitist selection: 0-no 1-yes\n"))
    type_mutation = int(input("Type mutation: 0-rare mutation 1-normal mutation\n"))

def calculate_lg():
    global lg
    dim = (B-A) * pow(10, precision)
    lg = int(math.log2(dim)) + 1

def from_base_2_to_10(individual):
    current_pow = 1
    number = 0
    global lg
    for i in range(lg - 1, -1, -1):
        number += individual[i] * current_pow
        current_pow *= 2
    number = (B-A)/(pow(2, lg) - 1) * number + A
    return number

#generate a new specimen with random genes
def generate_individual():
    individual = []
    for i in range(lg):
        gene = random.randint(0,1)
        individual.append(gene)
    return individual

def generate_population():
    global population
    g.write("Initial population:\n")
    for i in range(dim_population):
        individual = generate_individual()
        population.append(individual)
        g.write(f"{i + 1}: {individual}, x = {from_base_2_to_10(individual)} f ={calculate_fitness(individual)} \n")
    g.write("\n")

def calculate_fitness(individual):
    number = from_base_2_to_10(individual)
    value = a * number * number + b * number + c
    return value

def selection():

    global  intermediate_population
    intermediate_population = []
    fitness = []
    for i in range(len(population)):
        fitness.append(calculate_fitness(population[i]))  # save the fitnesses in an array because we will use it later
    total_fitness = sum(fitness)
    intervals = []
    maxim_fitness = 0

    global fittest_from_current_generation
    #global fittest_values
    fittest_from_current_generation = []
    global values_fittest

    #calculate the fittest element from current generation
    for i in range(len(population)):
        current_fitness = fitness[i]
        if current_fitness > maxim_fitness:
            maxim_fitness = current_fitness
            #fittest_from_current_generation = copy.deepcopy(population[i])
            fittest_from_current_generation = population[i].copy()

    values_fittest.append(calculate_fitness(fittest_from_current_generation)) # the vector with the fittest from all generation -->this will be ploted

    #calculate for each elem the probability to be selected current_fit/total_fit
    #generate the intervals
    previous = 0
    if ct == 1:
        g.write("Probabilities of selection:\n")
    for i in range(len(population)):
        current_probability = fitness[i] / total_fitness
        intervals.append(current_probability + previous)
        previous += current_probability
        if ct == 1:
            g.write(f"cromozom {i+1} probability {current_probability} \n")

    if ct == 1:
        g.write("\n Intervals probabilities of selection: \n")
        for interval in intervals:
            g.write(str(interval) + " ")
        g.write("\n")

    #dim_population times we will generate a random number and we will select the number that correspond to that interval
    for i in range(dim_population - is_elitist):
        value = random.random()
        index = binary_search(value, intervals)
        intermediate_population.append(population[index])
        if ct == 1:
            g.write(f"u = {value} select cromozom {index}\n")

    if ct == 1:
        g.write("After selection:\n")
        for i in range(len(intermediate_population)):
            g.write(f"{i+1}: {intermediate_population[i]}, x = {from_base_2_to_10(intermediate_population[i])} \n")

#for each element from intermediate_population generate a number and see if it will be selected for recombination
def crossing_over():
    if ct == 1:
        g.write(f"Probability of crossing-over: {p_cross}\n")

    global recombinated_population
    recombinated_population = []
    for_recombination = []

    for i in range(len(intermediate_population)):
        value = random.random()
        if(value <= p_cross):
            for_recombination.append(i)
            if ct == 1:
                g.write(f"{i+1}: {intermediate_population[i]} u = {value} <= {p_cross} participates\n")
        else:
            recombinated_population.append(intermediate_population[i])
            if ct == 1:
                g.write(f"{i + 1}: {intermediate_population[i]} u = {value}\n")

    # if just one -- do nothing -- just add it
    if len(for_recombination) == 1:
        recombinated_population.append(intermediate_population[for_recombination[0]])
    # recombinate circular for the odd case
    while(len(for_recombination) > 0 and len(for_recombination) != 1):
        if(len(for_recombination) == 3):
            index = random.randint(0, lg)#generate a random point breakage --> interchange until that point
            indx1 = for_recombination[0]
            indx2 = for_recombination[1]
            indx3 = for_recombination[2]

            elem1 = intermediate_population[indx1]
            elem2 = intermediate_population[indx2]
            elem3 = intermediate_population[indx3]

            first1 = elem1[:index]
            first2 = elem2[:index]
            first3 = elem3[:index]
            last1 = elem1[index:]
            last2 = elem2[index:]
            last3 = elem3[index:]

            if ct == 1:
                g.write(f"Recombination of cromozom {indx1} with cromozom {indx2} and {indx3}: \n {elem1} {elem2} {elem3}\n")

            elem1 = first1 + last3
            elem2 = first2 + last1
            elem3 = first3 + last2

            if ct == 1:
                g.write(f"Result: {elem1} {elem2} {elem3}\n")

            recombinated_population.append(elem1)
            recombinated_population.append(elem2)
            recombinated_population.append(elem3)
            for_recombination = []
        else:
            index = random.randint(0, lg) #to be sure that we will interchange smth --> can start with 1
            #choose last 2 elements
            #indx1 = for_recombination[len(for_recombination) - 1]
            indx1 = for_recombination.pop()
            elem1 = intermediate_population[indx1]

            #indx2 = for_recombination[len(for_recombination) - 1]
            indx2 = for_recombination.pop()
            elem2 = intermediate_population[indx2]

            first1 = elem1[:index]
            first2 = elem2[:index]

            last1 = elem1[index:]
            last2 = elem2[index:]

            if ct == 1:
                g.write(f"Recombination of cromozom {indx1} with cromozom {indx2}: \n {elem1} {elem2}\n")
            #concate and generate 2 new elements
            elem1 = first1 + last2
            elem2 = first2 + last1
            if ct == 1:
                g.write(f"Result: {elem1} {elem2}\n")

            recombinated_population.append(elem1)
            recombinated_population.append(elem2)

    if ct == 1:
        g.write("\nAfter recombination:\n")
        for i in range(len(recombinated_population)):
            g.write(f"{i + 1}. {recombinated_population[i]} x = {from_base_2_to_10(recombinated_population[i])}\n")

#for each gene in each cromozom we will generate a random number and if it's <= p_mutation --> we will change it
def normal_mutation():

    global mutated_population
    mutated_population = []
    mutated = set()
    for i in range(len(recombinated_population)):
        for j in range(lg):
            value = random.random()
            if value <= p_mutation:
                recombinated_population[i][j] = 1 - recombinated_population[i][j]
                mutated.add(i + 1)
        mutated_population.append(recombinated_population[i])

    # if elitist the fittest element is qualified
    if is_elitist:
        mutated_population.append(fittest_from_current_generation)

        if ct == 1:
            g.write(f"BEING ELITIST WE WILL SELECT AND THE FITTEST ELEM\n")

    if ct == 1:
        g.write(f"It have been modified cromozoms: \n")
        for cromozom in mutated:
            g.write(str(cromozom)+"\n")

        g.write("\nAfter mutation:\n")
        for i in range(len(mutated_population)):
            g.write(f"{i + 1}. {mutated_population[i]} x = {from_base_2_to_10(mutated_population[i])}\n")

def rare_mutation():
    global mutated_population
    mutated_population = []
    mutated = set()
    for i in range(len(recombinated_population)):
        value = random.random()
        if value <= p_mutation:
                j = random.randint(0, lg -1)
                recombinated_population[i][j] = 1 - recombinated_population[i][j]
                mutated.add(i + 1)
        mutated_population.append(recombinated_population[i])

    # if elitist the fittest element is qualified
    if is_elitist:
        mutated_population.append(fittest_from_current_generation)

        if ct == 1:
            g.write(f"BEING ELIIST WE WILL SELECT AND THE FITTEST ELEM\n")

    if ct == 1:
        g.write(f"It have been modified cromozoms: \n")
        for cromozom in mutated:
            g.write(str(cromozom) + "\n")

        g.write("\nAfter mutation:\n")
        for i in range(len(mutated_population)):
            g.write(f"{i + 1}. {mutated_population[i]} x = {from_base_2_to_10(mutated_population[i])}\n")

def binary_search(value, intervals):
    left = 0
    right = len(intervals) - 1
    while(left <= right):
        middle = (left + right) // 2
        if(middle == 0 and value < intervals[middle]):
            return 0
        if(value == intervals[middle] or (value < intervals[middle] and value > intervals[middle - 1])):
            return middle
        if (value < intervals[middle]):
            right = middle - 1
        else:
            left = middle + 1

def evolution():
    g.write("Maximum evolution: \n")
    for fit in values_fittest:
        g.write(str(fit) +"\n")
    plt.plot(values_fittest)
    plt.show()

def init_population():
    global population
    population = []

if __name__ == '__main__':

    g = open("output.txt", "w")
    read()
    calculate_lg()
    generate_population()
    ct = 1

    while ct <= steps:
        selection()
        crossing_over()

        if type_mutation == 0:
           normal_mutation()
        else:
            rare_mutation()

        init_population()
        for individual in mutated_population:
            population.append(individual)
        ct += 1
    evolution()
