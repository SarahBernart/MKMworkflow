import os
import shutil

# Define the temperature range
T_start = 0
T_end = 1000
T_step = 5

# Ensure temp folder exists
output_folder = "temp"
os.makedirs(output_folder, exist_ok=True)

# Path to ped.py script
ped_script = "ped.py"

# Read the ped.py script content
with open(ped_script, "r") as file:
    ped_content = file.readlines()

# Loop over the temperature range
for T in range(T_start, T_end + 1, T_step):
    # Modify temperature in ped.py
    modified_content = []
    for line in ped_content:
        if line.startswith("T ="):
            modified_content.append(f"T = {T}\n")
        else:
            modified_content.append(line)
    
    # Write modified ped.py
    modified_ped_path = "temp_ped.py"
    with open(modified_ped_path, "w") as file:
        file.writelines(modified_content)
    
    # Execute modified ped.py
    os.system(f"python {modified_ped_path}")
    
    # Rename and move ped.txt to temp folder
    output_file = os.path.join(output_folder, f"ped_{T}.txt")
    shutil.move("data/ped.txt", output_file)

print("All temperature calculations completed and saved in the 'temp' folder.")



