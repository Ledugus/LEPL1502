import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
from scipy.interpolate import interp1d


# Internal functions (not-used for user input)
def load_csv(file, delimiter=";"):
    # Load the CSV file
    df = pd.read_csv(
        file,
        delimiter=delimiter,
        skiprows=1,  # Adjust if the actual column names are in a different row
    )
    return df


def count_csv_columns(file_path, delimiter=";"):
    with open(file_path, "r") as file:
        # Read the first line of the file
        first_line = file.readline()
        # Count the number of columns based on the delimiter
        num_columns = (
            first_line.count(delimiter) + 1
        )  # Add 1 since there's one more column than delimiters
    return num_columns


# Graphing functions
def plot_from_csv(
    file_path,
    export_path,
    time_start,
    time_end,
    label_list,
    title,
    time_scale="$\mu s$",
    save_file=False,
    show_plot=True,
    return_data=False,
):
    """
    Plots data from a CSV file.

    Parameters:
    - file_path: str, the path to the CSV file to read.
    - export_path: str, the path where the plot should be saved if save_to_file is True.
    - time_start: float, the starting time value for the x-axis.
    - time_end: float, the ending time value for the x-axis.
    - labels: list of str, labels for each data series in the plot.
    - title: str, title of the plot
    - time_unit: str, the unit of time to display on the x-axis, default is microseconds.
    - save_to_file: bool, if True, the plot will be saved to the specified export_path.
    - show_plot: bool, if True, the plot will be displayed
    - return_data: bool, if True, returns data from csv file used to plot

    The function reads the specified CSV file, plots each column against the first column as time,
    and optionally saves the plot to a file. Assumes the first column is time and subsequent columns are data series.
    """
    # Handle errors
    if time_end < time_start:
        raise ValueError("time_end can't be smaller than time_start")

    # Open file
    df = load_csv(file_path, delimiter=",")

    # Count number of columns in file
    num_col = count_csv_columns(file_path, delimiter=",")
    # replace all infinite values with NaN
    df.replace(["∞", "-∞"], 1000, inplace=True)
    # Convert columns to appropriate data types, replace commas with dots, and convert to numeric if needed
    for i in range(num_col):
        df[df.columns[i]] = df[df.columns[i]].astype(float)
    # Plot the data
    plt.figure(figsize=(10, 6))

    # Plotting the original data
    try:
        for i in range(1, num_col):
            plt.plot(df[df.columns[0]], df[df.columns[i]], label=f"{label_list[i-1]}")
    # Handle error for inproper labeling
    except IndexError:
        raise IndexError(
            f"The number of labels for your plot should be '{num_col-1}', not '{len(label_list)}'."
        )

    plt.xlabel(f"Temps ({time_scale})", fontsize=18)
    plt.ylabel("Tension ($V$)", fontsize=18)
    plt.title(f"{title}", fontsize=24)
    plt.legend(fontsize=22, loc="upper right")

    plt.xlim(time_start, time_end)

    plt.tick_params(axis="both", labelsize=14)  # Increase the size of the tick labels

    # Save File
    if save_file is True:
        plt.savefig(
            export_path,
            format="pdf",
        )

    # Return data
    if return_data is True:
        return df

    # Show plot
    if show_plot is True:
        plt.show()


def plot_tangent(
    file_path,
    export_path,
    time_start,
    time_end,
    label_list,
    target_channel,
    time_scale="$\mu s$",
    tangent_point_index=0,
    dx=0.1,  # dx can now be non-integer
    tangent_width=0.1,
    save_file=False,
    show_plot=True,
):
    """
    Plots data from a CSV file and draws a tangent line at a specified point along a target data channel.

    The function first plots the data from the specified CSV file using plot_from_csv. It then
    interpolates the data for the target channel to calculate the slope of the tangent line at the point
    specified by tangent_point_index. A tangent line of specified width is then plotted over the data.

    Parameters:
    - file_path (str): The path to the CSV file containing the data to plot.
    - export_path (str): The path where the plot should be saved if save_file is True.
    - time_start (float): The starting time value for the x-axis of the plot.
    - time_end (float): The ending time value for the x-axis of the plot.
    - label_list (list of str): Labels for each data series in the plot.
    - target_channel (int): The index of the data column (channel) for which to draw the tangent line.
    - time_scale (str): The unit of time to display on the x-axis, default is microseconds.
    - tangent_point_index (int): The index of the point in the data at which to calculate and draw the tangent line.
    - dx (float): The delta x used for calculating the slope of the tangent, can be non-integer.
    - tangent_width (float): The width of the tangent line plotted on the graph.
    - save_file (bool): If True, the plot will be saved to the specified export_path.
    - show_plot (bool): If True, the plot will be displayed.

    Returns:
    None
    """
    df = plot_from_csv(
        file_path,
        export_path,
        time_start,
        time_end,
        label_list,
        time_scale,
        show_plot=False,
        return_data=True,
    )

    x = df.iloc[:, 0].values
    y = df.iloc[:, target_channel].values

    # Create an interpolation function
    interp_func = interp1d(x, y, kind="cubic")

    x0 = x[tangent_point_index]
    y0 = interp_func(x0)

    # Use the interpolation function to find y-values at x0 +/- dx for slope calculation
    x_plus_dx = x0 + dx
    y_plus_dx = interp_func(x_plus_dx)

    slope = (y_plus_dx - y0) / dx

    # Define the range for the tangent line, centered around x0
    x_line = np.linspace(x0, x0 + tangent_width, 100)

    # Equation of the tangent line
    tangent = slope * (x_line - x0) + y0

    plt.plot(x_line, tangent, color="red", label="Tangent Line")

    plt.tick_params(axis="both", which="major", labelsize=14)

    if show_plot:
        plt.legend()
        plt.show()

    if save_file:
        plt.savefig(export_path, format="pdf")


def plot_vertical_lines(
    x_value, y_value, transparency=0.5, color_x="Purple", color_y="Purple"
):
    plt.axvline(x=x_value, alpha=transparency, color=color_x, linestyle="--")
    plt.axhline(y=y_value, alpha=transparency, color=color_y, linestyle="--")
    plt.xticks([x_value])
    plt.yticks([y_value])


plot_from_csv(
    "data/V_R_f_18.27_R68_C104.csv",
    "edp.pdf",
    -50,
    50,
    [r"$V_1(t)$", r"$V_{R1}$"],
    "",
    time_scale="$\mu s$",
    save_file=False,
    show_plot=True,
)
