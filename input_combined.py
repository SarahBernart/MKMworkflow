import os

def extract_named_blocks_from_mkm(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()

    in_reaction_section = False
    current_label = None
    current_block = []
    named_blocks = {}

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("&reactions"):
            in_reaction_section = True
            continue
        elif stripped.startswith("&settings"):
            break

        if in_reaction_section:
            if stripped.startswith("#"):
                if current_block and current_label:
                    block_str = "\n".join(current_block).strip()
                    if block_str not in named_blocks.values():
                        named_blocks[current_label] = block_str
                    current_block = []
                current_label = stripped[1:].strip()
            elif stripped == "":
                if current_block and current_label:
                    block_str = "\n".join(current_block).strip()
                    if block_str not in named_blocks.values():
                        named_blocks[current_label] = block_str
                    current_block = []
                    current_label = None
            else:
                current_block.append(line.rstrip())

    if current_block and current_label:
        block_str = "\n".join(current_block).strip()
        if block_str not in named_blocks.values():
            named_blocks[current_label] = block_str

    return named_blocks

def extract_compounds_block(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()

    in_compound_section = False
    compound_header = []
    current_label = None
    current_block = []
    named_compound_blocks = {}

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("&compounds"):
            in_compound_section = True
            continue
        elif stripped.startswith("&reactions"):
            break

        if in_compound_section:
            if stripped.startswith("#"):
                if current_block and current_label:
                    block_str = "\n".join(current_block).strip()
                    if block_str not in named_compound_blocks.values():
                        named_compound_blocks[current_label] = block_str
                    current_block = []
                current_label = stripped[1:].strip()
            elif stripped == "":
                if current_block and current_label:
                    block_str = "\n".join(current_block).strip()
                    if block_str not in named_compound_blocks.values():
                        named_compound_blocks[current_label] = block_str
                    current_block = []
                    current_label = None
            elif current_label:
                current_block.append(line.rstrip())
            else:
                compound_header.append(line.rstrip())

    if current_block and current_label:
        block_str = "\n".join(current_block).strip()
        if block_str not in named_compound_blocks.values():
            named_compound_blocks[current_label] = block_str

    return compound_header, named_compound_blocks

def build_combined_input_with_compounds(mechanisms, output_dir="Kinetics_All", output_file="input.mkm"):
    compound_header = [
        "&compounds",
        "N2;   0;  0.9",
        "CO;   0;  0.08",
        "O2;   0;  0.02",
        "CO2;  0;  0",
        ""
    ]

    settings_block = """
&settings
PRESSURE = 1
TYPE = SEQUENCERUN
DEBUG = 1
NPAR=12
REAGENTS = {CO}, {O2}
KEYCOMPONENTS = {CO2}
DRC=1

&runs
400;1e6;1e-12;1e-12
450;1e6;1e-12;1e-12
500;1e6;1e-12;1e-12
550;1e6;1e-12;1e-12
600;1e6;1e-12;1e-12
650;1e6;1e-12;1e-12
700;1e6;1e-12;1e-12
750;1e6;1e-12;1e-12
800;1e6;1e-10;1e-10
850;1e6;1e-10;1e-10
900;1e6;1e-10;1e-10
950;1e6;1e-10;1e-10
1000;1e6;1e-10;1e-10
1050;1e6;1e-10;1e-10
1100;1e6;1e-10;1e-10
1150;1e6;1e-10;1e-10
1200;1e6;1e-10;1e-10
""".strip()

    all_compound_blocks = {}
    all_reaction_blocks = {}

    for mech in mechanisms:
        folder = f"Kinetics_{mech}"
        file_path = os.path.join(folder, "input.mkm")
        if os.path.isfile(file_path):
            _, compound_blocks = extract_compounds_block(file_path)
            for label, block in compound_blocks.items():
                if block not in all_compound_blocks.values():
                    all_compound_blocks[label] = block

            reaction_blocks = extract_named_blocks_from_mkm(file_path)
            for label, block in reaction_blocks.items():
                if block not in all_reaction_blocks.values():
                    all_reaction_blocks[label] = block
        else:
            print(f"⚠️ File not found: {file_path}")

    os.makedirs(output_dir, exist_ok=True)
    full_output_path = os.path.join(output_dir, output_file)

    with open(full_output_path, "w") as f:
        f.write("\n".join(compound_header) + "\n\n")
        for label, block in all_compound_blocks.items():
            f.write(f"# {label}\n")
            f.write(block + "\n\n")

        f.write("&reactions\n\n")
        for label, block in all_reaction_blocks.items():
            f.write(f"# {label}\n")
            f.write(block + "\n\n")

        f.write(settings_block + "\n")

    print(f"✅ Combined mechanism written to: {full_output_path}")

if __name__ == "__main__":
    mechanisms = ["L2", "Lh", "Mix", "Mvk", "Mix2"]
    build_combined_input_with_compounds(mechanisms)



