import itertools

macro_path = "/lustre/cmswork/coradin/GRID_SEARCH/200/"
param_range = [0.6, 0.85, 1] 

all_combinations = itertools.product(param_range, repeat=3)
i = 1
for combination in all_combinations:
    CFI0, CFI1, CF01 = combination
    file_content = f'''
TROOT::SetMacroPath("{macro_path}");
.L SNNT13.C
SNN_Tracking(100000, 10, 8, 10, "100k_100br_{i//5}.root", 1, {CFI0}, {CFI1}, {CF01});
.q
'''
    filename = f"starting_parameters_CF/start_{i}.cmd"
    i+=1
    with open(filename, "w") as file:
        file.write(file_content)
    print(f"File {filename} creato con successo.")

