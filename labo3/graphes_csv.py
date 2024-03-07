import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps


def graph(file, export_name, time_start, time_end, label1, label2, save_file=False):
    # Load the CSV file
    df = pd.read_csv(
        file,
        delimiter=";",
        skiprows=1,  # Adjust if the actual column names are in a different row
    )

    # Convert columns to appropriate data types, replace commas with dots, and convert to numeric if needed
    df[df.columns[0]] = df[df.columns[0]].str.replace(",", ".").astype(float)
    df[df.columns[1]] = df[df.columns[1]].str.replace(",", ".").astype(float)
    df[df.columns[2]] = df[df.columns[2]].str.replace(",", ".").astype(float)

    # Plot the data
    plt.figure(figsize=(10, 6))

    # Plotting the original data
    plt.plot(df[df.columns[0]], df[df.columns[1]], label=f"{label1}")
    plt.plot(df[df.columns[0]], df[df.columns[2]], label=f"{label2}")

    plt.xlabel("Temps (us)")
    plt.ylabel("Tension (V)")
    plt.title("Tensions en fonction du temps")
    plt.legend()

    plt.xlim(time_start, time_end)

    # Save File
    if save_file is True:
        plt.savefig(export_name, format="pdf")

    plt.show()


def graph_tan(
    file, export_name, time_start, time_end, slope, line_x_end, save_file=False
):
    # Load the CSV file
    df = pd.read_csv(
        file,
        delimiter=";",
        skiprows=1,  # Adjust if the actual column names are in a different row
    )

    # Convert columns to appropriate data types, replace commas with dots, and convert to numeric if needed
    df[df.columns[0]] = df[df.columns[0]].str.replace(",", ".").astype(float)
    df[df.columns[1]] = df[df.columns[1]].str.replace(",", ".").astype(float)
    df[df.columns[2]] = df[df.columns[2]].str.replace(",", ".").astype(float)

    # Plot the data
    plt.figure(figsize=(10, 6))

    # Plotting the original data
    plt.plot(df[df.columns[0]], df[df.columns[1]], label="$V(t)$")
    plt.plot(df[df.columns[0]], df[df.columns[2]], label="$V_c(t)$")

    # Plot horizontal limit line
    plt.axhline(y=1.75, linestyle="--", color="Green", alpha=0.5)

    # Plot vertical line to determine tau
    plt.axvline(x=20, linestyle="--", color="Blue", alpha=0.5)

    # Generate x-values for the line within a specific range
    line_x_values = np.linspace(0, line_x_end, 100)
    # Calculate y-values for the line using the slope
    line_y_values = slope * line_x_values

    # Plot the line with variable slope
    plt.plot(line_x_values, line_y_values, color="Red")

    # Get current ticks
    ticks = plt.xticks()[0]
    # Add the new tick while ensuring it's within the current x-limits
    if time_start < 20 < time_end:
        new_ticks = np.append(ticks, 20)
        plt.xticks(sorted(new_ticks))

    plt.xlabel("Temps (us)")
    plt.ylabel("Tension (V)")
    plt.title("Déterminer $\\tau$ graphiquement")
    plt.legend()

    plt.xlim(time_start, time_end)

    if save_file is True:
        plt.savefig(export_name, format="pdf")

    plt.show()


def graph_with_integration(
    file,
    export_name,
    time_start,
    time_end,
    integration_start,
    integration_end,
    save_file=False,
):
    # Load the CSV file
    df = pd.read_csv(
        file,
        delimiter=";",
        skiprows=1,  # Adjust if the actual column names are in a different row
    )

    # Convert columns to appropriate data types, replace commas with dots, and convert to numeric if needed
    df[df.columns[0]] = df[df.columns[0]].str.replace(",", ".").astype(float)
    df[df.columns[1]] = df[df.columns[1]].str.replace(",", ".").astype(float)
    df[df.columns[2]] = df[df.columns[2]].str.replace(",", ".").astype(float)

    # Plot the data
    plt.figure(figsize=(10, 6))

    # Plotting the original data
    plt.plot(df[df.columns[0]], df[df.columns[1]], label="$V(t)$")
    plt.plot(df[df.columns[0]], df[df.columns[2]], label="$V_c(t)$")

    # Shade the area under the curve for $V_c(t)$ within the integration period
    # Find indices corresponding to the integration limits
    integration_indices = (df[df.columns[0]] >= integration_start) & (
        df[df.columns[0]] <= integration_end
    )

    # Shade the area under $V_c(t)$
    plt.fill_between(
        df[df.columns[0]][integration_indices],
        0,
        df[df.columns[2]][integration_indices],
        color="lightgreen",
        alpha=0.5,
    )

    # Calculate the integral using Simpson's rule for $V_c(t)$
    integral_value = simps(
        df[df.columns[2]][integration_indices], df[df.columns[0]][integration_indices]
    )
    print(
        f"Integral value between {integration_start} and {integration_end} for $V_c(t)$: {integral_value}"
    )

    plt.xlabel("Temps (us)")
    plt.ylabel("Tension (V)")
    plt.title("Déterminer $V_cmoy(t)$")
    plt.legend()

    plt.xlim(time_start, time_end)

    if save_file is True:
        plt.savefig(export_name, format="pdf")

    plt.show()


graph(
    "labo3experience+-.csv",
    "exp2-2.pdf",
    -100,
    100,
    r"$V_{{in+}}$",
    r"$V_{{in-}}$",
    save_file=True,
)

graph(
    "labo3potentiometrea1,5.csv",
    "exp2-3.pdf",
    -100,
    100,
    r"$V_{{out}}$",
    r"$V_{{in-}}$",
    save_file=True,
)

# Example usage with integrating between x=10 and x=50
# graph_with_integration("S3100k100pf5khz_01.csv", "integration.pdf", -300, 300, 0, 200)
