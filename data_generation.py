import math
import os
import re
import shutil
import numpy as np

from nodal_info import (
    select_nodes_in_ellipsoid,
    damage_bdf,
    get_floating_nodes,
    get_box_nodes,
)


def _positive_normal(rng, mean, std, min_value=1e-6):
    value = rng.normal(mean, std)
    while value <= 0:
        value = rng.normal(mean, std)
    return max(value, min_value)


def _random_orientation(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def _bounded_normal(rng, mean, std, low, high):
    value = rng.normal(mean, std)
    while value < low or value > high:
        value = rng.normal(mean, std)
    return value


def _damage_number_from_input(rng, d_n):
    """
    Convert d_n into the actual number of damage locations for one sample.

    d_n can be either:
      - int, e.g. d_n=1 or d_n=3
      - range tuple/list, e.g. d_n=(1, 3), sampled inclusively from 1, 2, 3
    """
    if isinstance(d_n, (tuple, list, np.ndarray)):
        if len(d_n) != 2:
            raise ValueError("d_n range must have exactly two values, e.g. d_n=(1, 3).")
        low, high = int(d_n[0]), int(d_n[1])
        if low > high:
            raise ValueError("d_n range lower bound must be <= upper bound.")
        if low < 1:
            raise ValueError("d_n must be >= 1.")
        return int(rng.integers(low, high + 1))

    damage_number = int(d_n)
    if damage_number < 1:
        raise ValueError("d_n must be >= 1.")
    return damage_number


def write_multiple_damage_bdf(
    selected_nodes_groups,
    out_path,
    damage_ratio=0.5,
    base_bdf_path="base_file/FEM_only.bdf",
):
    """
    Write one BDF with multiple damaged regions using the same damage_ratio.

    This is intentionally separate from the existing damage_bdf/write_damage_bdf
    path.  For multiple regions, it creates only one damaged material/property
    pair per original property ID, instead of creating a new material for every
    damage location.

    selected_nodes_groups:
        Iterable of selected-node lists/sets.  Each item is one damage region.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    shutil.copy2(base_bdf_path, out_path)

    with open(out_path, "r", errors="ignore") as f:
        lines = f.readlines()

    def ff(line):
        return [line[i:i + 8].strip() for i in range(0, len(line.rstrip("\n")), 8)]

    def field(v):
        return str(v)[:8].rjust(8)

    def card(v):
        return str(v)[:8].ljust(8)

    def make(fields):
        return card(fields[0]) + "".join(field(x) for x in fields[1:]).rstrip() + "\n"

    def to_float(s):
        s = str(s).strip().replace("D", "E").replace("d", "e")
        if not s:
            return None
        if "e" in s.lower():
            return float(s)
        m = re.fullmatch(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)", s)
        if m:
            return float(m.group(1) + "e" + m.group(2))
        return float(s)

    def bdf_float_8(v):
        v = float(v)
        if not math.isfinite(v):
            raise ValueError(f"invalid BDF float: {v}")
        if v == 0:
            return "0."

        for prec in range(7, -1, -1):
            s = f"{v:.{prec}f}".rstrip("0").rstrip(".")
            if s in ("", "-0"):
                s = "0"
            if s == "0" and v != 0:
                continue
            if len(s) <= 8:
                return s

        exp = int(math.floor(math.log10(abs(v))))
        mant = v / (10 ** exp)
        for prec in range(6, -1, -1):
            m = f"{mant:.{prec}f}".rstrip("0").rstrip(".")
            s = f"{m}{exp:+d}"
            if len(s) <= 8:
                return s

        raise ValueError(f"could not format {v} as an 8-character BDF float")

    node_groups = []
    for selected_nodes in selected_nodes_groups:
        group = set(map(int, selected_nodes))
        if group:
            node_groups.append(group)

    if not node_groups:
        return out_path, []

    elem_cards = {
        "CBAR": 2,
        "CBEAM": 2,
        "CROD": 2,
        "CTRIA3": 3,
        "CQUAD4": 4,
        "CTETRA": 4,
        "CHEXA": 8,
    }
    prop_cards = {"PSHELL", "PBAR"}

    props = {}
    mats = {}
    selected_elems = []

    for line in lines:
        flds = ff(line)
        c = flds[0]

        if c in elem_cards:
            eid = int(flds[1])
            pid = int(flds[2])
            conn = [int(x) for x in flds[3:3 + elem_cards[c]] if x]

            # Important difference from simply using one unioned node set:
            # an element must be fully inside at least one damage region.
            if conn and any(all(x in group for x in conn) for group in node_groups):
                selected_elems.append((eid, pid))

        elif c in prop_cards:
            props[int(flds[1])] = {"fields": flds.copy(), "mid": int(flds[2])}

        elif c == "MAT1":
            mats[int(flds[1])] = flds.copy()

    if not selected_elems:
        return out_path, []

    next_pid = max(props.keys()) + 1000
    next_mid = max(mats.keys()) + 1000

    pid_map = {}
    add_lines = []

    # Only one damaged material/property pair is created per original property.
    # This avoids extra duplicate materials when several locations share one damage_ratio.
    for old_pid in sorted(set(pid for _, pid in selected_elems)):
        if old_pid not in props:
            continue
        old_mid = props[old_pid]["mid"]
        if old_mid not in mats:
            continue

        new_pid = next_pid
        next_pid += 1
        new_mid = next_mid
        next_mid += 1

        mf = mats[old_mid].copy()
        while len(mf) < 7:
            mf.append("")
        mf[1] = str(new_mid)

        e_val = to_float(mf[2])
        if e_val is not None:
            mf[2] = bdf_float_8(e_val * damage_ratio)

        g_val = to_float(mf[3])
        if g_val is not None:
            mf[3] = bdf_float_8(g_val * damage_ratio)

        pf = props[old_pid]["fields"].copy()
        while len(pf) < 8:
            pf.append("")
        pf[1] = str(new_pid)
        pf[2] = str(new_mid)

        pid_map[old_pid] = new_pid
        add_lines.append("$ Damaged material/property\n")
        add_lines.append(make(mf))
        add_lines.append(make(pf))

    selected_eid_set = set(eid for eid, _ in selected_elems)

    out_lines = []
    inserted = False

    for line in lines:
        flds = ff(line)
        c = flds[0]

        if not inserted and c == "GRID":
            out_lines.extend(add_lines)
            inserted = True

        if c in elem_cards:
            eid = int(flds[1])
            if eid in selected_eid_set:
                old_pid = int(flds[2])
                if old_pid in pid_map:
                    line = line[:16] + field(pid_map[old_pid]) + line[24:]

        out_lines.append(line)

    if not inserted:
        out_lines.extend(add_lines)

    with open(out_path, "w") as f:
        f.writelines(out_lines)

    print("selected elements:", len(selected_elems))
    print("new properties:", len(pid_map))
    print("pid map:", pid_map)
    print("saved:", out_path)

    return out_path, [eid for eid, _ in selected_elems]


def generate_data_sets(n=10, base_dir="base_file", temp_dir="Temp", info_dir="_data", seed=None, d_n=1):
    rng = np.random.default_rng(seed)
    os.makedirs(temp_dir, exist_ok=True)
    floating_nodes = get_floating_nodes(info_dir)
    box_nodes = get_box_nodes(info_dir)
    excluded_nodes = floating_nodes | box_nodes
    all_nodes = np.loadtxt(os.path.join(info_dir, "nodes.csv"), delimiter=",", skiprows=1, usecols=[0], dtype=int)
    candidate_nodes = [node_id for node_id in all_nodes.tolist() if node_id not in excluded_nodes]

    if not candidate_nodes:
        raise ValueError("No candidate nodes found after excluding box and floating nodes.")

    for i in range(1, n + 1):
        subfolder = os.path.join(temp_dir, f"sample_{i}")
        os.makedirs(subfolder, exist_ok=True)

        for name in os.listdir(base_dir):
            src = os.path.join(base_dir, name)
            dst = os.path.join(subfolder, name)
            if os.path.isfile(src):
                shutil.copy2(src, dst)

        damage_number = _damage_number_from_input(rng, d_n)
        if damage_number > len(candidate_nodes):
            raise ValueError(
                f"damage_number={damage_number} is larger than available candidate nodes "
                f"({len(candidate_nodes)})."
            )

        # Keep the single-damage case as close as possible to the old random path.
        if damage_number == 1:
            nodes = [int(rng.choice(candidate_nodes))]
        else:
            # Choose damage centers without replacement so one sample has n distinct locations.
            nodes = [int(x) for x in rng.choice(candidate_nodes, size=damage_number, replace=False)]

        damage_specs = []
        for node in nodes:
            a = _positive_normal(rng, 30, 5)
            b = _positive_normal(rng, 30, 5)
            c = _positive_normal(rng, 30, 5)
            ome = _random_orientation(rng)
            selected_nodes = select_nodes_in_ellipsoid(
                node=node,
                a=a,
                b=b,
                c=c,
                Ome=ome,
                data_dir=info_dir,
            )
            damage_specs.append(
                {
                    "node": node,
                    "a": a,
                    "b": b,
                    "c": c,
                    "Ome": ome,
                    "selected_nodes": selected_nodes,
                }
            )

        damage_ratio = np.clip(rng.normal(0.5, 0.1), 0.1, 0.7)

        with open(os.path.join(subfolder, "params.txt"), "w") as f:
            # Keep the old keys for the single-damage case, so older post-processing
            # scripts that read damaged_node/a/b/c/Ome/damage_ratio still work.
            if damage_number == 1:
                spec = damage_specs[0]
                f.write(f"damaged_node: {spec['node']}\n")
                f.write(f"a: {spec['a']:.12g}\n")
                f.write(f"b: {spec['b']:.12g}\n")
                f.write(f"c: {spec['c']:.12g}\n")
                f.write("Ome:\n")
                f.write(np.array2string(spec["Ome"], separator=", "))
                f.write("\n")
                f.write(f"damage_ratio: {damage_ratio:.12g}\n")
                f.write("damage_number: 1\n")
            else:
                f.write(f"damage_number: {damage_number}\n")
                f.write("damaged_nodes: " + ", ".join(str(s["node"]) for s in damage_specs) + "\n")
                f.write(f"damage_ratio: {damage_ratio:.12g}\n")

                for j, spec in enumerate(damage_specs, start=1):
                    f.write(f"damage_{j}_damaged_node: {spec['node']}\n")
                    f.write(f"damage_{j}_a: {spec['a']:.12g}\n")
                    f.write(f"damage_{j}_b: {spec['b']:.12g}\n")
                    f.write(f"damage_{j}_c: {spec['c']:.12g}\n")
                    f.write(f"damage_{j}_Ome:\n")
                    f.write(np.array2string(spec["Ome"], separator=", "))
                    f.write("\n")
                    f.write(f"damage_{j}_selected_nodes_count: {len(spec['selected_nodes'])}\n")

        out_bdf_path = os.path.join(subfolder, "FEM_only.bdf")
        base_bdf_path = os.path.join(base_dir, "FEM_only.bdf")

        if damage_number == 1:
            damage_bdf(
                selected_nodes=damage_specs[0]["selected_nodes"],
                out_path=out_bdf_path,
                damage_ratio=damage_ratio,
                base_bdf_path=base_bdf_path,
            )
        else:
            write_multiple_damage_bdf(
                selected_nodes_groups=[spec["selected_nodes"] for spec in damage_specs],
                out_path=out_bdf_path,
                damage_ratio=damage_ratio,
                base_bdf_path=base_bdf_path,
            )


if __name__ == "__main__":
    generate_data_sets(n=10)
