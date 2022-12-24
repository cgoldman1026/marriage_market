from sympy import * # import the required package
from sympy.abc import x, y, z, a, b
init_printing(use_unicode=True) # make the printing look nice

preference_matrix = [
    [0.3, 0.5],
    [-0.2, -0.2]
]

demographic_matrix = [
    [100, 100],
    [100, 100]
]

male_types = set(["r", "p"])
female_types = set(["r", "p"])

all_equations = []

# for every type of male in our male types
for male_idx, male_type in enumerate(male_types):
    # for every type of female in our female types
    for female_idx, female_type in enumerate(female_types):
        # print debugging info
        print(f"male_idx = {male_idx}, male_type = {male_type}")
        print(f"female_idx = {female_idx}, female_type = {female_type}")

        # make the preference symbols
        u1, u2, u3 = symbols(f"u_{male_type}_{female_type} u_{male_type}_0 u_0_{female_type}")

        # get the alpha value
        alpha_value = preference_matrix[male_idx][female_idx]

        # make the equation components
        numerator = 2 * u1
        denominator = sqrt(u2 * u3)
        result = exp(alpha_value)

        # put the equation together
        final_eq = Eq(numerator / denominator, result)

        # add equation to the list of equations
        all_equations.append(final_eq)
        


# for every type of male in our male types
for male_idx, male_type in enumerate(male_types):

    # collect a list of symbols
    symbol_list = []

    # make the first male symbol
    symbol_list.append(symbols(f"u_{male_type}_0"))

    # for every type of female in our female types
    for female_idx, female_type in enumerate(female_types):
        # make a female symbol
        symbol_list.append(symbols(f"u_{male_type}_{female_type}"))
    
    # build the male constraint equation
    final_eq = Eq(demographic_matrix[male_idx][female_idx] - sum(symbol_list[1:]), symbol_list[0])
    
    # add equation to the list of equations
    all_equations.append(final_eq)


# for every type of female in our female types
for female_idx, female_type in enumerate(female_types):
    
    # collect a list of symbols
    symbol_list = []

    # make the first female symbol
    symbol_list.append(symbols(f"u_0_{female_type}"))

    # for every type of male in our male types
    for male_idx, male_type in enumerate(male_types):
        # make a female symbol
        symbol_list.append(symbols(f"u_{male_type}_{female_type}"))
    
    # build the female constraint equation
    final_eq = Eq(demographic_matrix[male_idx][female_idx] - sum(symbol_list[1:]), symbol_list[0])
    
    # add equation to the list of equations
    all_equations.append(final_eq)



for equation in all_equations:
    print(equation, "\n")
        
        

solutions = solve(all_equations, dict=True)
print(solutions)

# x, y, z = symbols("x y z")

# formulas = [
#     Eq(x**2 + y**2 + z**2, 6),
#     Eq(x**2 - y**2 + 2*z**2, 2),
#     Eq(2*x**2 + y**2 - z**2, 3)
# ]
# solutions = solve(formulas, dict=True)
# for solution in solutions:
#     print(solution)


