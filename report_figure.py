import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection


def _require_cv2():
    try:
        import cv2
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "OpenCV is required for drawing report circles. "
            "Install opencv-python or call plot_element_strain_presence only."
        ) from exc

    return cv2


def draw_dashed_circle(
    img,
    center,
    radius,
    color=(0, 0, 0),
    thickness=4,
    dash_length=18,
    gap_length=12,
):
    """
    Draw a dashed circle using short arc segments.

    center: (x, y)
    radius: circle radius in pixels
    color: BGR color, e.g. black = (0, 0, 0)
    """
    cv2 = _require_cv2()

    circumference = 2 * np.pi * radius
    dash_angle = 360 * dash_length / circumference
    gap_angle = 360 * gap_length / circumference

    angle = 0
    while angle < 360:
        start_angle = angle
        end_angle = min(angle + dash_angle, 360)

        cv2.ellipse(
            img,
            center,
            (radius, radius),
            angle=0,
            startAngle=start_angle,
            endAngle=end_angle,
            color=color,
            thickness=thickness,
            lineType=cv2.LINE_AA,
        )

        angle += dash_angle + gap_angle


def add_circles_to_image(input_path, output_path):
    cv2 = _require_cv2()
    img = cv2.imread(str(input_path))

    if img is None:
        raise FileNotFoundError(f"Could not read image: {input_path}")

    # ------------------------------------------------------------
    # Circle positions
    # Coordinate system:
    # x increases to the right
    # y increases downward
    #
    # You can manually change these later.
    # These exact same coordinates are used for both figures.
    # ------------------------------------------------------------

    solid_circle_center = (180, 520)     # left-middle black circle
    solid_circle_radius = 90

    dashed_circle_1_center = (470, 455)
    dashed_circle_1_radius = 95

    dashed_circle_2_center = (770, 600)
    dashed_circle_2_radius = 105

    # Solid black circle
    cv2.circle(
        img,
        solid_circle_center,
        solid_circle_radius,
        color=(0, 0, 0),
        thickness=5,
        lineType=cv2.LINE_AA,
    )

    # Dashed black circles
    draw_dashed_circle(
        img,
        dashed_circle_1_center,
        dashed_circle_1_radius,
        color=(0, 0, 0),
        thickness=5,
        dash_length=20,
        gap_length=14,
    )

    draw_dashed_circle(
        img,
        dashed_circle_2_center,
        dashed_circle_2_radius,
        color=(0, 0, 0),
        thickness=5,
        dash_length=20,
        gap_length=14,
    )

    cv2.imwrite(str(output_path), img)


def draw_report_circles(image_paths=None, output_paths=None):
    """Draw the report annotation circles on the configured images."""
    if image_paths is None:
        image_paths = [
            Path("img/undamaged.png"),
            Path("img/70.png"),
        ]

    if output_paths is None:
        output_paths = [
            Path("img/undamaged_with_circles.png"),
            Path("img/70_with_circles.png"),
        ]

    if len(image_paths) != len(output_paths):
        raise ValueError("image_paths and output_paths must have the same length.")

    for input_path, output_path in zip(image_paths, output_paths):
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        add_circles_to_image(input_path, output_path)

    return [Path(path) for path in output_paths]


def _set_axes_equal(ax):
    """Force equal axis scale in a 3D Matplotlib axis."""
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    max_range = max(x_range, y_range, z_range) / 2.0

    x_middle = sum(x_limits) / 2.0
    y_middle = sum(y_limits) / 2.0
    z_middle = sum(z_limits) / 2.0

    ax.set_xlim3d([x_middle - max_range, x_middle + max_range])
    ax.set_ylim3d([y_middle - max_range, y_middle + max_range])
    ax.set_zlim3d([z_middle - max_range, z_middle + max_range])


def _resolve_strain_column(strain_df, strain_column):
    strain_columns = [
        column for column in strain_df.columns
        if str(column).strip().upper() != "ID"
    ]

    if not strain_columns:
        raise ValueError("No strain columns found after the ID column.")

    if isinstance(strain_column, int):
        try:
            return strain_columns[strain_column]
        except IndexError as exc:
            raise ValueError(
                f"strain_column index {strain_column} is out of range. "
                f"Available strain columns: 0..{len(strain_columns) - 1}"
            ) from exc

    if strain_column not in strain_df.columns:
        raise ValueError(f"Column not found in strain CSV: {strain_column}")

    return strain_column


def plot_element_strain_presence(
    data_dir="_data",
    strain_path="Temp/sample_1/1.csv",
    strain_column=0,
    output_path="img/strain_presence.png",
    zero_tol=0.0,
    nonzero_color="#d62728",
    zero_color="#b8b8b8",
    show=False,
):
    """
    Plot elements with nonzero strain in one color and zero/blank strain in another.

    strain_column is a zero-based index into the strain columns, excluding ID.
    The default 0 means the first strain value column after ID.
    """
    data_dir = Path(data_dir)
    strain_path = Path(strain_path)

    nodes_df = pd.read_csv(data_dir / "nodes.csv")
    elements_df = pd.read_csv(data_dir / "elements.csv")
    strain_df = pd.read_csv(strain_path, low_memory=False)

    selected_column = _resolve_strain_column(strain_df, strain_column)

    valid_ids = pd.to_numeric(strain_df["ID"], errors="coerce")
    strain_values = pd.to_numeric(strain_df[selected_column], errors="coerce")

    valid_strain_df = pd.DataFrame({
        "elem_id": valid_ids,
        "strain": strain_values,
    }).dropna(subset=["elem_id"])
    valid_strain_df["elem_id"] = valid_strain_df["elem_id"].astype(int)

    elements_with_strain = elements_df.merge(
        valid_strain_df,
        on="elem_id",
        how="left",
    )
    raw_strain = elements_with_strain["strain"]
    elements_with_strain["plot_strain"] = raw_strain.fillna(0.0)
    elements_with_strain["has_nonzero_strain"] = (
        elements_with_strain["plot_strain"].abs() > zero_tol
    )

    node_xyz = {
        int(row["node_id"]): (float(row["x"]), float(row["y"]), float(row["z"]))
        for _, row in nodes_df.iterrows()
    }
    node_cols = [column for column in elements_df.columns if column.startswith("n")]

    line_segments = {
        True: [],
        False: [],
    }
    face_segments = {
        True: [],
        False: [],
    }

    for _, elem in elements_with_strain.iterrows():
        conn = []
        for column in node_cols:
            value = elem[column]
            if pd.isna(value):
                continue
            try:
                node_id = int(value)
            except ValueError:
                continue
            if node_id in node_xyz:
                conn.append(node_id)

        if len(conn) < 2:
            continue

        points = [node_xyz[node_id] for node_id in conn]
        has_nonzero = bool(elem["has_nonzero_strain"])
        elem_type = str(elem["type"]).strip().upper()

        if elem_type in {"CTRIA3", "CQUAD4"} and len(points) >= 3:
            face_segments[has_nonzero].append(points)
        elif elem_type in {"CBAR", "CBEAM", "CROD"}:
            line_segments[has_nonzero].append(points[:2])
        else:
            line_segments[has_nonzero].append(points + [points[0]])

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")

    if face_segments[False]:
        zero_faces = Poly3DCollection(
            face_segments[False],
            facecolor=zero_color,
            edgecolor=zero_color,
            linewidth=0.08,
            alpha=0.22,
        )
        ax.add_collection3d(zero_faces)

    if line_segments[False]:
        zero_lines = Line3DCollection(
            line_segments[False],
            colors=zero_color,
            linewidths=0.25,
            alpha=0.25,
        )
        ax.add_collection3d(zero_lines)

    if face_segments[True]:
        nonzero_faces = Poly3DCollection(
            face_segments[True],
            facecolor=nonzero_color,
            edgecolor=nonzero_color,
            linewidth=0.08,
            alpha=0.85,
        )
        ax.add_collection3d(nonzero_faces)

    if line_segments[True]:
        nonzero_lines = Line3DCollection(
            line_segments[True],
            colors=nonzero_color,
            linewidths=0.35,
            alpha=0.85,
        )
        ax.add_collection3d(nonzero_lines)

    all_points = nodes_df[["x", "y", "z"]].to_numpy(float)
    ax.auto_scale_xyz(all_points[:, 0], all_points[:, 1], all_points[:, 2])
    _set_axes_equal(ax)

    nonzero_count = int(elements_with_strain["has_nonzero_strain"].sum())
    zero_count = int(len(elements_with_strain) - nonzero_count)
    blank_count = int(raw_strain.isna().sum())
    numeric_zero_count = int(
        (
            raw_strain.notna()
            & (raw_strain.abs() <= zero_tol)
        ).sum()
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(
        "Element Strain Presence\n"
        f"{selected_column} | nonzero={nonzero_count}, zero/blank={zero_count}"
    )
    ax.legend(
        handles=[
            Patch(color=nonzero_color, label="Nonzero strain"),
            Patch(color=zero_color, label="Zero or blank strain"),
            Line2D([], [], color="none", label=f"Blank values: {blank_count}"),
        ],
        loc="upper right",
    )

    plt.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300)

    if show:
        plt.show()

    summary = {
        "selected_column": selected_column,
        "total_elements": int(len(elements_with_strain)),
        "nonzero_elements": nonzero_count,
        "zero_or_blank_elements": zero_count,
        "blank_elements": blank_count,
        "numeric_zero_elements": numeric_zero_count,
        "nonzero_by_type": elements_with_strain.loc[
            elements_with_strain["has_nonzero_strain"], "type"
        ].value_counts().to_dict(),
        "zero_or_blank_by_type": elements_with_strain.loc[
            ~elements_with_strain["has_nonzero_strain"], "type"
        ].value_counts().to_dict(),
        "output_path": str(output_path) if output_path is not None else None,
    }

    return fig, ax, summary


if __name__ == "__main__":
    _, _, strain_summary = plot_element_strain_presence()

    print(f"Saved strain presence plot: {strain_summary['output_path']}")
