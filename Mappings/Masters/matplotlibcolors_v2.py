import matplotlib.colors as mcolors
from matplotlib import cm


def redefine_colors():
    # Redefine color shortcuts
    mcolors.ColorConverter.colors['b'] = "#1171BE"  # Change 'b' (blue) to a new color
    mcolors.ColorConverter.colors['o'] = "#DD5400"  # Change 'r' (red) to a new color
    mcolors.ColorConverter.colors['g'] = "#3BAA32"  # Change 'g' (green) to a new color
    mcolors.ColorConverter.colors['y'] = "#EDB120"  # Change 'y' (yellow) to a new color
    mcolors.ColorConverter.colors['p'] = "#7E2F8E"  # Change 'm' (magenta) to a new color
    mcolors.ColorConverter.colors['c'] = "#2FBEEF"  # Change 'c' (cyan) to a new color
    mcolors.ColorConverter.colors['r'] = "#A2142F"  # Change 'r' (red) to a new color
    mcolors.ColorConverter.colors['m'] = "#D1038B"  # Change 'r' (red) to a new color

if __name__ != "__main__":
    redefine_colors()