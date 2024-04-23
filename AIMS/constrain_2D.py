with open("geometry.in", "r") as f:
    lines = [line.strip() for line in f.readlines()]

idx = 0
dim = 0
while True:
    if lines[idx].strip().startswith("lattice"):
        lattice = lines[idx].split()[1:]
        if dim in [0, 1]:
            lines.insert(idx+1, "constrain_relaxation z")
        elif dim == 2:
            lines.insert(idx+1, "constrain_relaxation .true.")
        else:
            break
    idx += 1
    if idx == len(lines):
        break


with open("geometry.in", "w") as f:
    f.write("\n".join(lines))