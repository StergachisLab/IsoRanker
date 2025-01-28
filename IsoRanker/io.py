import pandas as pd

def load_data(filename):
    """Load isoform data from a CSV file."""
    return pd.read_csv(filename, sep=",")
