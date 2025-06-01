import sys

def parse_msh(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    nodes = {}
    elements = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # ----------- Parse $Nodes Section -----------
        if line == "$Nodes":
            i += 1
            header = list(map(int, lines[i].strip().split()))
            num_entity_blocks, total_num_nodes = header[0], header[1]
            i += 1

            for _ in range(num_entity_blocks):
                entity_header = list(map(int, lines[i].strip().split()))
                if len(entity_header) != 4:
                    raise ValueError(f"Unexpected node block header: {lines[i]}")
                entity_dim, entity_tag, parametric, num_nodes_in_block = entity_header
                i += 1

                node_ids = []
                while len(node_ids) < num_nodes_in_block:
                    node_ids.extend(map(int, lines[i].strip().split()))
                    i += 1

                for nid in node_ids:
                    coords = list(map(float, lines[i].strip().split()))
                    if len(coords) != 3:
                        raise ValueError(f"Invalid coordinate line: {lines[i]}")
                    nodes[nid] = coords
                    i += 1
            continue

        # ----------- Parse $Elements Section -----------
        if line == "$Elements":
            i += 1
            header = list(map(int, lines[i].strip().split()))
            num_entity_blocks, total_num_elements = header[0], header[1]
            i += 1

            for _ in range(num_entity_blocks):
                entity_header = list(map(int, lines[i].strip().split()))
                if len(entity_header) != 4:
                    raise ValueError(f"Unexpected element block header: {lines[i]}")
                entity_dim, entity_tag, element_type, num_elements_in_block = entity_header
                i += 1

                for _ in range(num_elements_in_block):
                    tokens = list(map(int, lines[i].strip().split()))
                    element_id = tokens[0]
                    conn = tokens[1:]

                    if element_type == 2 and len(conn) == 3:  # 2 = triangle
                        elements.append(conn)

                    i += 1
            continue

        i += 1

    return nodes, elements

def renumber_nodes(nodes):
    sorted_ids = sorted(nodes.keys())
    id_map = {old_id: new_id for new_id, old_id in enumerate(sorted_ids)}
    coords = [nodes[old_id] for old_id in sorted_ids]
    return id_map, coords

def write_petsc(filename, coords, elements, id_map):
    with open(filename, 'w') as f:
        f.write(f"{len(elements)} {len(coords)}\n")
        for elem in elements:
            mapped = [id_map[nid] for nid in elem]
            f.write(" ".join(map(str, mapped)) + "\n")
        for pt in coords:
            f.write(f"{pt[0]} {pt[1]} {pt[2]}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 gmsh2petsc.py input.msh output.dat")
        return

    msh_file = sys.argv[1]
    dat_file = sys.argv[2]

    try:
        nodes, elements = parse_msh(msh_file)
    except Exception as e:
        print(f"Error while parsing {msh_file}: {e}")
        return

    id_map, coords = renumber_nodes(nodes)
    write_petsc(dat_file, coords, elements, id_map)
    print(f"✅ Converted {msh_file} to {dat_file} successfully.")

if __name__ == "__main__":
    main()
