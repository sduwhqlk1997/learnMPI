#!/usr/bin/env python3
"""
Convert Gmsh-4.x ASCII .msh to two PETSc-style files:

  <prefix>.is   # 单元三元组 | 节点类型 | 边界边二元组
  <prefix>.vec  # 节点坐标 (x y z)

节点类型：0 = interior, 1 = Neumann, 2 = Dirichlet
依赖：纯标准库，兼容 Python ≥3.8
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Dict

# ---------------------------------------------------------------------------#
# 1. 解析 $PhysicalNames
# ---------------------------------------------------------------------------#
def read_physical(lines: List[str], idx: int) -> Tuple[Dict[int, str], int]:
    n = int(lines[idx + 1])
    tag2name: Dict[int, str] = {}
    for j in range(n):
        dim, tag, name = lines[idx + 2 + j].split(maxsplit=2)
        tag2name[int(tag)] = name.strip('"')
    return tag2name, idx + 2 + n + 1        # +1 跳过 $EndPhysicalNames

# ---------------------------------------------------------------------------#
# 2. 解析 $Entities —— 记录 “几何曲线 tag → 物理 tag”
# ---------------------------------------------------------------------------#
def read_entities(lines: List[str], idx: int) -> Tuple[Dict[int, int], int]:
    pt, crv, srf, vol = map(int, lines[idx + 1].split())
    i = idx + 2 + pt                       # 跳过 point-entities
    curve2phys: Dict[int, int] = {}
    for _ in range(crv):
        nums = list(map(float, lines[i].split()))
        tag = int(nums[0])
        n_phys = int(nums[7])
        first_phys = int(nums[8]) if n_phys else 0
        curve2phys[tag] = first_phys
        i += 1
    # 跳过 surface / volume 描述
    i += srf + vol
    return curve2phys, i + 1               # +1 跳过 $EndEntities

# ---------------------------------------------------------------------------#
# 3. 解析 $Nodes —— 支持“多行 ID + 多行坐标 + numNodes=0”
# ---------------------------------------------------------------------------#
def read_nodes(lines: List[str], idx: int):
    num_blocks, total_nodes, *_ = map(int, lines[idx + 1].split())
    i = idx + 2
    coords: Dict[int, Tuple[float, float, float]] = {}
    for _ in range(num_blocks):
        ent_dim, ent_tag, param, nnod = map(int, lines[i].split())
        i += 1
        # 读取 node IDs
        node_ids: List[int] = []
        while len(node_ids) < nnod:
            node_ids.extend(map(int, lines[i].split()))
            i += 1
        # 读取坐标 (parametric==0 ⇒ 每节点 3 个浮点数)
        xyz: List[float] = []
        need = (3 if param == 0 else 6) * nnod
        while len(xyz) < need:
            xyz.extend(map(float, lines[i].split()))
            i += 1
        for k, nid in enumerate(node_ids):
            coords[nid] = tuple(xyz[3 * k : 3 * k + 3])
    return coords, i + 1                   # +1 跳过 $EndNodes

# ---------------------------------------------------------------------------#
# 4. 解析 $Elements —— 单元 / 边 + 节点类型
# ---------------------------------------------------------------------------#
def read_elements(
    lines: List[str],
    idx: int,
    curve2phys: Dict[int, int],
    tag2name: Dict[int, str],
):
    etype_nnodes = {1: 2, 2: 3}            # 1=line, 2=triangle
    num_blocks, *_ = map(int, lines[idx + 1].split())
    i = idx + 2
    triangles: List[Tuple[int, int, int]] = []
    boundary_edges: List[Tuple[int, int]] = []
    node_type = defaultdict(int)           # default 0 (interior)

    for _ in range(num_blocks):
        ent_dim, ent_tag, etype, nelem = map(int, lines[i].split())
        i += 1
        nnodes = etype_nnodes.get(etype, 0)
        for _ in range(nelem):
            data = list(map(int, lines[i].split()))
            node_tags = data[-nnodes:]
            if etype == 1:                 # line → 边界
                boundary_edges.append(tuple(node_tags))
                ptag = curve2phys.get(ent_tag, 0)
                name = tag2name.get(ptag, "")
                for n in node_tags:
                    if name == "dirichlet":
                        node_type[n] = 2
                    elif name == "neumann" and node_type[n] != 2:
                        node_type[n] = 1
            elif etype == 2:               # triangle
                triangles.append(tuple(node_tags))
            i += 1
    return triangles, boundary_edges, node_type, i + 1  # +1 跳过 $EndElements

# ---------------------------------------------------------------------------#
# 5. 主解析函数
# ---------------------------------------------------------------------------#
def parse_msh(file: Path):
    lines = file.read_text().splitlines()
    i = 0
    tag2name, curve2phys = {}, {}
    coords = {}
    triangles = []
    boundary_edges = []
    node_type = defaultdict(int)

    while i < len(lines):
        s = lines[i].strip()
        if s == "$PhysicalNames":
            tag2name, i = read_physical(lines, i)
        elif s == "$Entities":
            curve2phys, i = read_entities(lines, i)
        elif s == "$Nodes":
            coords, i = read_nodes(lines, i)
        elif s == "$Elements":
            triangles, boundary_edges, node_type, i = read_elements(
                lines, i, curve2phys, tag2name
            )
        else:
            i += 1

    # 重新编号：Gmsh → 连续局部索引
    gid = sorted(coords)
    lid = {g: l for l, g in enumerate(gid)}
    tri_local  = [(lid[a], lid[b], lid[c]) for a, b, c in triangles]
    edge_local = [(lid[a], lid[b]) for a, b in boundary_edges]
    types      = [node_type[g] for g in gid]
    xyz        = [coords[g] for g in gid]
    return tri_local, types, edge_local, xyz

# ---------------------------------------------------------------------------#
# 6. 写输出
# ---------------------------------------------------------------------------#
def write_output(prefix: str, cells, types, edges, xyz):
    with open(prefix + ".is", "w") as f:
        for a, b, c in cells:
            f.write(f"{a} {b} {c}\n")
        f.write("\n")
        for t in types:
            f.write(f"{t}\n")
        f.write("\n")
        for a, b in edges:
            f.write(f"{a} {b}\n")

    with open(prefix + ".vec", "w") as f:
        for x, y, _ in xyz:
            f.write(f"{x} {y}\n")

# ---------------------------------------------------------------------------#
# 7. CLI
# ---------------------------------------------------------------------------#
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mshfile", help="Gmsh 4.x ASCII .msh file")
    ap.add_argument("-o", "--out", default="mesh", help="output file prefix")
    ns = ap.parse_args()

    tri, tp, edg, crd = parse_msh(Path(ns.mshfile))
    write_output(ns.out, tri, tp, edg, crd)
    print(f"✔ 生成 {ns.out}.is  和  {ns.out}.vec")

if __name__ == "__main__":
    main()
