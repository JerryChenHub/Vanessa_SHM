import argparse
import os
import shutil
import json
import hashlib
from pathlib import Path
import numpy as np

from .paths import BASE_DIR, DATA_DIR, CASES_DIR, require_working_output
from .nodal_info import (
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


def write_multiple_damage_bdf(selected_nodes_groups, out_path, damage_ratio=0.3, base_bdf_path=BASE_DIR / "FEM_only.bdf"):
    from .bdf_damage import write_damage_regions
    return write_damage_regions(selected_nodes_groups, out_path, damage_ratio, base_bdf_path)


def generate_data_sets(n=10, base_dir=BASE_DIR, temp_dir=CASES_DIR, info_dir=DATA_DIR, seed=None, d_n=1):
    """Copy base inputs into workspace cases and apply stiffness damage.

    ``temp_dir`` remains the output-directory argument for existing Python callers.
    """
    temp_dir = require_working_output(temp_dir)
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
        if any(Path(subfolder).glob('sol144.f06')):
            raise FileExistsError(f"Refusing to overwrite a solved case: {subfolder}")
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
            _, modified_elements = damage_bdf(
                selected_nodes=damage_specs[0]["selected_nodes"],
                out_path=out_bdf_path,
                damage_ratio=damage_ratio,
                base_bdf_path=base_bdf_path,
            )
        else:
            _, modified_elements = write_multiple_damage_bdf(
                selected_nodes_groups=[spec["selected_nodes"] for spec in damage_specs],
                out_path=out_bdf_path,
                damage_ratio=damage_ratio,
                base_bdf_path=base_bdf_path,
            )

        if not modified_elements:
            raise ValueError(f"No elements changed in {subfolder}; do not label this as damaged.")
        manifest = {
            "schema_version": 1,
            "sample_index": i,
            "seed": seed,
            "stiffness_scale": float(damage_ratio),
            "damage_definition": "Scale MAT1 E/G for all active PSHELL membrane, bending, transverse shear and coupling material references, and for PBAR beams; preserve density and geometry.",
            "selection_rule": "All element nodes must be inside the same ellipsoid.",
            "modified_element_ids": sorted(modified_elements),
            "modified_element_count": len(modified_elements),
            "regions": [{"center_node": s["node"], "axes": [s["a"],s["b"],s["c"]],
                         "orientation": s["Ome"].tolist(),
                         "selected_node_ids": sorted(map(int,s["selected_nodes"]))} for s in damage_specs],
            "base_input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in Path(base_dir).iterdir() if p.suffix.lower() in {".bdf", ".inc"}},
            "case_input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in Path(subfolder).iterdir() if p.suffix.lower() in {".bdf", ".inc"}},
        }
        Path(subfolder, "case_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Generate damage cases from core/base into workspace/cases.")
    parser.add_argument("--samples", type=int, default=10)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--damage-regions", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=CASES_DIR)
    parser.add_argument("--base-dir", type=Path, default=BASE_DIR)
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    args = parser.parse_args()
    if args.samples < 1 or args.damage_regions < 1:
        parser.error("--samples and --damage-regions must be positive")
    generate_data_sets(args.samples, args.base_dir, args.output_dir, args.data_dir, args.seed, args.damage_regions)


if __name__ == "__main__":
    main()
