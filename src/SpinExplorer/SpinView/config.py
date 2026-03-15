import sys
from enum import Enum, auto

if sys.platform == "linux":
    platform = "linux"
    height = 30
elif sys.platform == "darwin":
    platform = "mac"
    height = 20
else:
    platform = "windows"
    height = 30

colours = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]

colour_options = [
            "Blue",
            "Orange",
            "Green",
            "Red",
            "Purple",
            "Brown",
            "Pink",
            "Gray",
            "Lime",
            "Turquoise         ",
            "Black",
            "Navy",
            "Tan",
            "Light Coral",
            "Maroon",
            "Light Green",
            "Deep Pink",
            "Fuchsia",
        ]

twoD_colours = [
            "#e41a1c",
            "#377eb8",
            "#4daf4a",
            "#984ea3",
            "#ff7f00",
            "#ff33eb",
            "tan",
            "lightcoral",
            "maroon",
            "lightgreen",
            "deeppink",
            "fuchsia",
        ]

reference_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "5.0",
            "10.0",
            "50.0",
        ]

multiply_range_values = [
            "1.01",
            "1.1",
            "1.5",
            "2",
            "5",
            "10",
            "50",
            "100",
            "1000",
            "10000",
            "100000",
            "1000000",
            "10000000",
            "100000000",
            "1000000000",
        ]

vertical_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "10.0",
            "50.0",
            "100.0",
            "1000.0",
            "10000",
        ]


class ScrollMode(Enum):
    ZOOM = auto()
    CONTOUR = auto()
    PLANE = auto()

def cycle_scroll_mode(current_mode,available_modes):
    current_index = available_modes.index(current_mode)
    return available_modes[(current_index+1)%len(available_modes)]

ONED_SCROLL_MODES = [ScrollMode.ZOOM]
TWOD_SCROLL_MODES = [ScrollMode.ZOOM, ScrollMode.CONTOUR]
THREED_SCROLL_MODES = [ScrollMode.ZOOM, ScrollMode.CONTOUR, ScrollMode.PLANE]
