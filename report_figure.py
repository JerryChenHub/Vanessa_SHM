import cv2
import numpy as np
from pathlib import Path


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


# Original figures
image_paths = [
    Path("img/undamaged.png"),
    Path("70.png"),
]

# Output figures
output_paths = [
    Path("img/undamaged_with_circles.png"),
    Path("70_with_circles.png"),
]

for input_path, output_path in zip(image_paths, output_paths):
    add_circles_to_image(input_path, output_path)

print("Finished drawing circles.")