import numpy as np
from scipy.optimize import fsolve, least_squares
import math
import time


def read_data(filename):

    with open(f"data/{filename}.csv", "r") as infile:
        # load the entire spreadsheet into an array
        spreadsheet = infile.readlines()

        # number of male types
        all_male_types = set([int(row.split(",")[0]) for row in spreadsheet[1:]]) - {0}
        num_male_types = len(all_male_types)

        # number of female types
        all_female_types = set([int(row.split(",")[2]) for row in spreadsheet[1:]]) - {0}
        num_female_types = len(all_female_types)

    return spreadsheet, num_male_types, num_female_types


def build_equation_templates():

    # print("vars = ", vars)
    
    # initialize a variable to hold all equations
    equation_templates = []
    
    # initialize a list to hold all initial guesses
    initial_guesses = []

    # initialize a list to store the bounds for each row
    lower_bounds = []
    upper_bounds = []


    # for each enumerated row in the spreadsheet
    for row_idx, row_data in enumerate(spreadsheet):

        # print("row_idx = ", row_idx)
        # print("row_data = ", row_data)

        # if this is the header row
        if row_idx == 0:
            continue

        # split the data of the row
        male_id, male_num, female_id, female_num, initial_guess, alpha = row_data.strip().split(",")
        
        # add the initial guess to the list of guesses
        initial_guesses.append(float(initial_guess))

        # find the bound for this row
        try:
            male_num = int(male_num)
        except:
            male_num = float("inf")
        try:
            female_num = int(female_num)
        except:
            female_num = float("inf")
        
        lower_bounds.append(0)
        upper_bounds.append(min(male_num, female_num))

        # turn the ids into integers so we can operate on them
        male_id = int(male_id)
        female_id = int(female_id)

        # if this is a normal equation row
        if male_id != 0 and female_id != 0:
            # make a note of everything we need for this equation
            new_equation_dict = {
                "numerator_idx": row_idx - 1,
                "denominator_idx_1": num_male_types * num_female_types + male_id - 1,
                "denominator_idx_2": num_male_types * num_female_types + num_male_types + female_id - 1,
                "exponential_piece": math.exp(float(alpha))
            }

            # add the equation to our list
            equation_templates.append(new_equation_dict)

        # if this is a male constraint equation
        elif female_id == 0:
            # make a note of everything we need for this equation
            new_equation_dict = {
                "starting_sum_idx": (male_id - 1) * num_female_types + 1,
                "ending_sum_idx": (male_id - 1) * num_female_types + num_female_types,
                "male_num": int(male_num),
                "current_row_idx": row_idx - 1
            }

            # add the equation to our list
            equation_templates.append(new_equation_dict)

        # if this is a female constraint equation
        elif male_id == 0:

            # make a note of everything we need for this equation
            new_equation_dict = {
                "starting_sum_idx": 1,
                "ending_sum_idx": num_male_types,
                "female_id": female_id,
                "num_female_types": num_female_types,
                "female_num": int(female_num),
                "current_row_idx": row_idx - 1
            }

            # add the equation to our list
            equation_templates.append(new_equation_dict)

    return equation_templates, initial_guesses, lower_bounds, upper_bounds



def func(vars):
    
    
    # for the index of each equation we need to build
    for equation_idx in range(len(all_equations)):
        
        # if this is a normal equation row
        if equation_idx < num_male_types * num_female_types:

            # print("normal equation with index = ", equation_idx)
            # for k in equation_templates[equation_idx]:
            #     print(f"{k} = {equation_templates[equation_idx][k]}")

            # build the numerator based on the dictionary
            numerator = 2 * vars[equation_templates[equation_idx]["numerator_idx"]]
            # build the denominator based on the dictionary
            denominator = math.sqrt((vars[equation_templates[equation_idx]["denominator_idx_1"]]) * (vars[equation_templates[equation_idx]["denominator_idx_2"]]))
            # update the new equation
            all_equations[equation_idx] = numerator / denominator - equation_templates[equation_idx]["exponential_piece"]

        # if this is a male constraint equation
        elif equation_idx < num_male_types * num_female_types + num_male_types:

            # accumulate the sum of the variables
            var_sum = 0

            # for each variable from starting to ending index
            for var_idx in range(equation_templates[equation_idx]["starting_sum_idx"], equation_templates[equation_idx]["ending_sum_idx"] + 1):
                # add the variable to the sum
                var_sum += vars[var_idx - 1]
            
            # update the new equation
            all_equations[equation_idx] = equation_templates[equation_idx]["male_num"] - var_sum - vars[equation_templates[equation_idx]["current_row_idx"]]

        # if this is a female constraint equation
        else:
            
            # accumulate the sum of the variables
            var_sum = 0

            # for each variable from starting index to ending index
            for var_idx in range(equation_templates[equation_idx]["starting_sum_idx"], equation_templates[equation_idx]["ending_sum_idx"] + 1):
                # add the variable to the sum
                var_sum += vars[(equation_templates[equation_idx]["female_id"]) + (var_idx - 1) * (equation_templates[equation_idx]["num_female_types"]) - 1]

            # update the new equation
            all_equations[equation_idx] = equation_templates[equation_idx]["female_num"] - var_sum - vars[equation_templates[equation_idx]["current_row_idx"]]
    

    # print("vars = ", vars)

    return all_equations



# name of the file you want to read data from

# all_names = []

infile = "example_data_frac500"
print("infile = ", infile)

time1 = time.time()

# read the data, get the number of each type
spreadsheet, num_male_types, num_female_types = read_data(infile)

time2 = time.time()

# build the equation templates and get the initial guesses from the spreadsheet
equation_templates, vars, lower_bounds, upper_bounds = build_equation_templates()

time3 = time.time()

# start with blank equations
all_equations = [0.0 for i in range((num_male_types * num_female_types) + num_male_types + num_female_types)]

time4 = time.time()

# # run the solver
# root = least_squares(func, vars, bounds=[lower_bounds, upper_bounds], xtol=1e-3)
root = least_squares(func, vars, method="trf", bounds=[lower_bounds, upper_bounds], xtol=1e-3)


time5 = time.time()

with open(f"solutions/{infile}.csv", "w") as outfile:
    existing_headers = spreadsheet[0].strip()
    existing_data = spreadsheet[1:]

    new_headers = f"{existing_headers},mu\n"
    outfile.write(new_headers)

    for idx in range(len(existing_data)):
        existing_row = existing_data[idx].strip()
        new_row = f"{existing_row},{str(root.x[idx])}\n"
        outfile.write(new_row)


read_data_time = time2 - time1
build_template_time = time3 - time2
full_prep_time = time4 - time1
solver_time = time5 - time4
overall_time = time5 - time1
solver_time_frac = (solver_time / overall_time) * 100
print("read_data_time = ", read_data_time)
print("build_template_time = ", build_template_time)
print("full_prep_time = ", full_prep_time)
print("solver_time = ", solver_time)
print("overall_time = ", overall_time)
print("solver_time_frac = ", solver_time_frac)








