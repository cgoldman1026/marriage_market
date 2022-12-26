import numpy as np
from scipy.optimize import fsolve, least_squares
import math
import time

outfile = open("julia_equations.jl", "w")

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
        
        # # add the initial guess to the list of guesses
        # initial_guesses.append(float(initial_guess))

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

            outfile.write(f"""\tF[{row_idx}] = (x[{row_idx}]) / (sqrt(x[{num_male_types * num_female_types + male_id}] * x[{num_male_types * num_female_types + num_male_types + female_id}])) -  exp({float(alpha)})\n""")

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

            outstring = f"""\tF[{row_idx}] = {int(male_num)}"""

            for var_idx in range((male_id - 1) * num_female_types + 1, (male_id - 1) * num_female_types + num_female_types + 1):
                outstring += f" - x[{var_idx}]"

            outstring += f" - x[{row_idx}]"

            outfile.write(outstring + "\n")

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

            outstring = f"""\tF[{row_idx}] = {int(female_num)}"""

            for var_idx in range(1, num_male_types + 1): 
                outstring += f""" - x[{(female_id) + (var_idx - 1) * (num_female_types)}]"""
            
            outstring += f""" - x[{row_idx}]"""

            outfile.write(outstring + "\n")
    


    return equation_templates, initial_guesses, lower_bounds, upper_bounds


infile = "example_data_frac40_guess"
print("infile = ", infile)

time1 = time.time()

# read the data, get the number of each type
spreadsheet, num_male_types, num_female_types = read_data(infile)

time2 = time.time()

# build the equation templates and get the initial guesses from the spreadsheet
equation_templates, vars, lower_bounds, upper_bounds = build_equation_templates()

outfile.close()





