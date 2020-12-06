import re
from pathlib import Path

def main():
    instance_list = Path("data/instance_list.txt")
    with instance_list.open() as f:
        data = [line[:-1] for line in f.readlines()]
    
    data_root = Path("data")
    for tsp in data:
        file = data_root / tsp
        with file.open() as f:
            tsp_data = f.readlines()
        
        dim = int(re.search(r"DIMENSION.+?(\d+)", "".join(tsp_data))[1])
        for i, line in enumerate(tsp_data):
            if "NODE_COORD_SECTION"  in line:
                break

        tsp_data.insert(i, f"MIN_NODE_NUM : {int(dim*0.75)}\n")

        with file.open("w") as f:
            f.writelines(tsp_data)
        

if __name__ == "__main__":
    main()