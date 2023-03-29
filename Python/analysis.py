import numpy as np
import matplotlib.pyplot as plt


def process_experiment1():
    # Load the CNOT counts from the Gaussian elimination and PMH algorithms.
    results_gauss = np.loadtxt(
        './resources/scalability_gauss.csv',
        dtype=int, delimiter=','
    )
    results_pmh = np.loadtxt(
        './resources/scalability_pmh.csv',
        dtype=int, delimiter=','
    )

    # Number of qubits and circuit size for each set of CNOT counts.
    qubit_range = results_gauss[:, 0]
    size_range = results_gauss[:, 1]

    # Get the raw CNOT counts from the two algorithms.
    counts_gauss = results_gauss[:, 2:]
    counts_pmh = results_pmh[:, 2:]

    # Calculate the mean circuit sizes.
    mean_gauss = np.mean(counts_gauss, axis=1)
    mean_pmh = np.mean(counts_pmh, axis=1)

    # Error bars are given by the minimum and maximum circuit sizes.
    bars_gauss = np.row_stack([
        mean_gauss - np.min(counts_gauss, axis=1),
        np.max(counts_gauss, axis=1) - mean_gauss
    ])
    bars_pmh = np.row_stack([
        mean_pmh - np.min(counts_pmh, axis=1),
        np.max(counts_pmh, axis=1) - mean_pmh
    ])

    # Create a new figure.
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111)

    # Plot the number of random row operations used to generate the parity matrices.
    ax.plot(qubit_range, size_range, color='black', linestyle='dashed')

    # Plot the mean circuit size versus number of qubits for the Gaussian elimination algorithm.
    ax.errorbar(
        qubit_range, mean_gauss, bars_gauss,
        label='Gauss', marker='o', ms=3, elinewidth=1, capsize=1.5
    )

    # Plot the mean circuit size versus number of qubits for the PMH algorithm.
    ax.errorbar(
        qubit_range, mean_pmh, bars_pmh,
        label='PMH', marker='o', ms=3, elinewidth=1, capsize=1.5
    )

    # Add x- and y-axis labels and create a legend.
    ax.set_xlabel('Number of Qubits', weight='bold', fontsize=10)
    ax.set_ylabel('Mean Circuit Size', weight='bold', fontsize=10)
    ax.grid(True)
    ax.legend()

    # Save the figure to the resources directory.
    fig.savefig('./resources/scalability.png', bbox_inches='tight', dpi=300)


def process_experiment2():
    # Range of qubit values considered.
    n_min = 10
    n_max = 230
    n_step = 20

    # Range of section sizes considered.
    m_min = 2
    m_max = 8
    m_step = 1

    # Load the raw CNOT counts from experiment 2.
    results = np.loadtxt(
        './resources/optimisation_raw.csv',
        dtype=int, delimiter=','
    )

    # Create a table of all zeros.
    table = np.zeros((
        ((n_max - n_min) // n_step) + 1,
        ((m_max - m_min) // m_step) + 1
    ))

    # Map n and m values to table indicies.
    results[:, 0] = (results[:, 0] - n_min) / n_step
    results[:, 1] = (results[:, 1] - m_min) / m_step

    # Calculate the mean circuit size for each (n, m) tuple.
    mean_counts = np.mean(results[:, 2:], axis=1)

    # Add the mean circuit sizes to the table.
    for i, (qubit, section_size) in enumerate(results[:, 0:2]):
        table[qubit, section_size] = mean_counts[i]

    # Add a column on the left for the values of n.
    table = np.column_stack((
        np.arange(n_min, n_max + n_step, n_step),
        table
    ))

    # Add a row at the top for the values of m.
    table = np.row_stack((
        np.arange(m_min - m_step, m_max + m_step, m_step),
        table
    ))

    # Set a dummy value for the top left-most cell.
    table[0, 0] = -1

    # Save the formatted results to a CSV file.
    np.savetxt(
        './resources/optimisation_processsed.csv', table,
        delimiter=',', fmt='%.1f'
    )


if __name__ == '__main__':
    # Process the raw data from experiments 1 and 2 (see the C code).
    process_experiment1()
    process_experiment2()
