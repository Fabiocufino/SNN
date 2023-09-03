import itertools

macro_path = "/lustre/cmswork/coradin/GRID_SEARCH/200/"
param_range = range(6, 11)  # Valori compresi tra 6 e 10

all_combinations = itertools.product(param_range, repeat=2)
i = 1
for combination in all_combinations:
    param_value_3, param_value_4 = combination
    file_content = f'''
TROOT::SetMacroPath("{macro_path}");
.L SNNT13.C
SNN_Tracking(100000, 10, {param_value_3}, {param_value_4}, "100k_200br_{i/5}.root", 1);
.q
'''
    filename = f"starting_parameters_neurons/start_{i}.cmd"
    i+=1
    with open(filename, "w") as file:
        file.write(file_content)
    print(f"File {filename} creato con successo.")

