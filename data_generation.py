import os
import shutil
import numpy as np

from nodal_info import (
    select_nodes_in_ellipsoid,
    damage_bdf,
    get_floating_nodes,
    get_box_nodes,
)


def _positive_normal(rng, mean, std, min_value=1e-6):
    # Sample a positive normal value.
    value = rng.normal(mean, std)
    while value <= 0:
        value = rng.normal(mean, std)
    return max(value, min_value)


def generate_data_sets(n=10, base_dir="base_file", temp_dir="Temp", info_dir="_data", seed=None):
    # Generate damaged FEM samples in Temp.
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

        node = int(rng.choice(candidate_nodes))
        radius = _positive_normal(rng, 30, 5)
        decay_rate = _positive_normal(rng, 0.5, 0.1)

        selected_nodes = select_nodes_in_ellipsoid(
            node=node,
            a=radius,
            b=radius,
            c=radius,
            Ome=np.eye(3),
            decay_rate=decay_rate,
            data_dir=info_dir,
        )

        damage_bdf(
            selected_nodes=selected_nodes,
            out_path=os.path.join(subfolder, "FEM_only.bdf"),
            damage_ratio=0.1,
            base_bdf_path=os.path.join(base_dir, "FEM_only.bdf"),
        )


if __name__ == "__main__":
    generate_data_sets(n=10)
