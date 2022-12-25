import numpy as np
from scipy.optimize import fsolve, least_squares
import math

infile_name = "example_data-2"
with open(f"{infile_name}.csv", "r") as infile:
    # load the entire spreadsheet into an array
    spreadsheet = infile.readlines()

    # number of male types
    all_male_types = set([int(row.split(",")[0]) for row in spreadsheet[1:]]) - {0}
    num_male_types = len(all_male_types)
    # # print("all_male_types = ", all_male_types)

    # number of female types
    all_female_types = set([int(row.split(",")[2]) for row in spreadsheet[1:]]) - {0}
    num_female_types = len(all_female_types)
    # # print("all_female_types = ", all_female_types)


def build_equations(vars):
    # print("vars = ", vars)
    # make a list to accumulate all equations into
    all_equations = []

    # for each enumerated row in the spreadsheet
    for row_idx, row_data in enumerate(spreadsheet):

        # print("row_idx = ", row_idx)
        # print("row_data = ", row_data)

        # if this is the header row
        if row_idx == 0:
            continue

        # split the data of the row
        male_id, male_num, female_id, female_num, alpha = row_data.strip().split(",")

        # turn the ids into integers so we can operate on them
        male_id = int(male_id)
        female_id = int(female_id)
        

        # if this is a normal equation row
        if male_id != 0 and female_id != 0:
            # print("Normal Equation Row: ")
            # turn the alpha number into a float so we can use it in the equation
            alpha = float(alpha)
            # print("alpha = ", alpha)

            # build the numerator piece
            numerator = 2 * vars[row_idx - 1]
            # print("row_idx - 1 = ", row_idx - 1)
            
            # build the denominator piece
            # print("(vars[num_male_types * num_female_types + male_id - 1]) = ", (vars[num_male_types * num_female_types + male_id - 1]))
            # print("(vars[num_male_types * num_female_types + num_male_types + female_id - 1]) = ", (vars[num_male_types * num_female_types + num_male_types + female_id - 1]))
            denominator = math.sqrt((vars[num_male_types * num_female_types + male_id - 1]) * (vars[num_male_types * num_female_types + num_male_types + female_id - 1]))
            
            
            
            # print("num_male_types * num_female_types + male_id - 1 = ", num_male_types * num_female_types + male_id - 1)
            # print("vars[num_male_types * num_female_types + male_id - 1] = ", vars[num_male_types * num_female_types + male_id - 1])
            # print("num_male_types * num_female_types + num_male_types + female_id - 1 = ", num_male_types * num_female_types + num_male_types + female_id - 1)
            # print("vars[num_male_types * num_female_types + num_male_types + female_id - 1] = ", vars[num_male_types * num_female_types + num_male_types + female_id - 1], "\n\n")
            # print("num_male_types * num_female_types + male_id - 1 = ", num_male_types * num_female_types + male_id - 1)
            # print("num_male_types * num_female_types + num_male_types + female_id - 1 = ", num_male_types * num_female_types + num_male_types + female_id - 1)
            
            # build the exponential piece
            exponential_piece = math.exp(alpha)
            # print("exponential_piece = ", exponential_piece)
            
            # put all the pieces together to get the new equation
            new_equation = numerator / denominator - exponential_piece
            all_equations.append(new_equation)

        # if this is a male constraint equation
        elif female_id == 0:
            # print("Male Constraint Equation Row: ")

            # get the starting index of the sum to take
            starting_sum_idx = (male_id - 1) * num_female_types + 1
            # print("starting_sum_idx = ", starting_sum_idx)

            # get the ending index of the sum to take
            ending_sum_idx = (male_id - 1) * num_female_types + num_female_types
            # print("ending_sum_idx = ", ending_sum_idx)

            # accumulate the sum of the variables
            var_sum = 0

            # for each variable from starting index to ending index
            for var_idx in range(starting_sum_idx, ending_sum_idx + 1):
                # print("var_idx = ", var_idx)
                
                # add the variable to the sum
                var_sum += vars[var_idx - 1]

            # turn the male_num into an integer
            male_num = int(male_num)
            # print("male_num = ", male_num)

            # build the equation
            new_equation = male_num - var_sum - vars[row_idx - 1]
            # print("row_idx - 1 = ", row_idx - 1)
            all_equations.append(new_equation)

        # if this is a female constraint equation
        elif male_id == 0:
            # # print("female constraint equation line")

            # get the starting index of the sum to take
            starting_sum_idx = 1
            # print("starting_sum_idx = ", starting_sum_idx)

            # get the ending index of the sum to take
            ending_sum_idx = num_male_types
            # print("ending_sum_idx = ", ending_sum_idx)

            # accumulate the sum of the variables
            var_sum = 0

            # for each variable from starting index to ending index
            for var_idx in range(starting_sum_idx, ending_sum_idx + 1):
                # add the variable to the sum
                var_sum += vars[female_id + (var_idx - 1) * num_female_types - 1]
                # print("female_id + (var_idx - 1) * num_female_types - 1 = ", female_id + (var_idx - 1) * num_female_types - 1)
            
            # turn the female number into an integer
            female_num = int(female_num)
            # print("female_num = ", female_num)

            # build the equation
            new_equation = female_num - var_sum - vars[row_idx - 1]
            # print("row_idx - 1 = ", row_idx - 1)
            all_equations.append(new_equation)
        
        # print("\n\n\n")

    return all_equations


vars = []
bds = []
for i in range((num_male_types * num_female_types)):
    vars.append(0)
    bds.append((0, 100))

for i in range(num_male_types + num_female_types):
    vars.append(50)
    bds.append((0, 100))

# print(bds)
# create an array that is going to hold all of our variables
# vars = [0.00001 for i in range((num_male_types * num_female_types) + num_male_types + num_female_types)]
# build_equations(vars)
# solve the system
root = least_squares(build_equations, vars, bounds = ((0, 100)))

with open(f"{infile_name}_solution.csv", "w") as outfile:
    existing_headers = spreadsheet[0].strip()
    existing_data = spreadsheet[1:]

    new_headers = f"{existing_headers},mu\n"
    outfile.write(new_headers)

    for idx in range(len(existing_data)):
        existing_row = existing_data[idx].strip()
        new_row = f"{existing_row},{str(root.x[idx])}\n"
        outfile.write(new_row)









