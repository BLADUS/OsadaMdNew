def extract_coordinates(input_file, output_file):
    with open(input_file, 'r') as f_input, open(output_file, 'w') as f_output:
        for line in f_input:
            if line.startswith("r1"):
                coordinates = line.split("=")[-1].strip().strip("()").split("; ")
                if len(coordinates) >= 2:
                    x, y = coordinates[0], coordinates[1]
                    f_output.write(f"r1:{x}; {y}\n")
            elif line.startswith("r2"):
                coordinates = line.split("=")[-1].strip().strip("()").split("; ")
                if len(coordinates) >= 2:
                    x, y = coordinates[0], coordinates[1]
                    f_output.write(f"r2:{x}; {y}\n")

# Пример использования:
input_file = "Osada_MD_7.txt"
output_file = "pars_output.txt"
extract_coordinates(input_file, output_file)
