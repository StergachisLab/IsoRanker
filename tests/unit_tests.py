import unittest
import pandas as pd
import numpy as np
from IsoRanker import (
    load_data,
    filter_based_on_counts,
    calculate_ranks_for_sample,
    NMD_test_statistic,
    Noncyclo_Expression_Outlier_LOE,
    Noncyclo_Expression_Outlier_GOE,
    Cyclo_Expression_Outlier_LOE,
    Cyclo_Expression_Outlier_GOE,
    NMD_rare_steady_state_transcript,
    calculate_z_score,
    apply_hypothesis_test,
    process_hypothesis_test,
    create_expression_matrix,
    create_long_format
)

# Private import helper functions
from IsoRanker.expression_matrix import parse_read_stats

class TestIO(unittest.TestCase):
    def test_load_data(self):
        test_csv = "test_data.csv"
        pd.DataFrame({"A": [1, 2], "B": [3, 4]}).to_csv(test_csv, index=False)
        df = load_data(test_csv)
        self.assertTrue(isinstance(df, pd.DataFrame))
        self.assertEqual(df.shape, (2, 2))

class TestPreprocessing(unittest.TestCase):
    def test_filter_based_on_counts(self):
        data = pd.DataFrame({
            "Isoform": ["iso1", "iso1", "iso2"],
            "cyclo_count": [5, 15, 20],
            "noncyclo_count": [25, 10, 5]
        })
        filtered_df = filter_based_on_counts(data, count_threshold=10)
        self.assertEqual(len(filtered_df), 3)

class TestRanking(unittest.TestCase):
    def test_calculate_ranks_for_sample(self):
        data = pd.DataFrame({
            "Sample": ["S1", "S1", "S2"],
            "z_score": [1.5, 2.5, 3.5],
            "test_statistic": [10, 20, 30],
            "Isoform": ["iso1", "iso2", "iso3"]
        })
        ranked_df = calculate_ranks_for_sample(data, "Isoform")
        self.assertIn("rank_top_99_5_percentile", ranked_df.columns)

class TestTestStatistic(unittest.TestCase):
    def test_NMD_test_statistic(self):
        data = pd.DataFrame({
            "Cyclo_TPM": [10, 20],
            "Noncyclo_TPM": [5, 15],
            "cyclo_count": [100, 200],
            "noncyclo_count": [50, 150]
        })
        result = NMD_test_statistic(data)
        self.assertIn("test_statistic", result.columns)

    def test_Noncyclo_Expression_Outlier_LOE(self):
        data = pd.DataFrame({
            "Noncyclo_TPM": [5, 15, 25],
            "cyclo_count": [10, 20, 30],
            "noncyclo_count": [5, 15, 25]
        })
        result = Noncyclo_Expression_Outlier_LOE(data)
        self.assertIn("test_statistic", result.columns)

    def test_Noncyclo_Expression_Outlier_GOE(self):
        data = pd.DataFrame({
            "Noncyclo_TPM": [5, 15, 25],
            "cyclo_count": [10, 20, 30],
            "noncyclo_count": [5, 15, 25]
        })
        result = Noncyclo_Expression_Outlier_GOE(data)
        self.assertIn("test_statistic", result.columns)

    def test_Cyclo_Expression_Outlier_LOE(self):
        data = pd.DataFrame({
            "Cyclo_TPM": [5, 15, 25],
            "cyclo_count": [10, 20, 30],
            "noncyclo_count": [5, 15, 25]
        })
        result = Cyclo_Expression_Outlier_LOE(data)
        self.assertIn("test_statistic", result.columns)

    def test_Cyclo_Expression_Outlier_GOE(self):
        data = pd.DataFrame({
            "Cyclo_TPM": [5, 15, 25],
            "cyclo_count": [10, 20, 30],
            "noncyclo_count": [5, 15, 25]
        })
        result = Cyclo_Expression_Outlier_GOE(data)
        self.assertIn("test_statistic", result.columns)

    def test_NMD_rare_steady_state_transcript(self):
        data = pd.DataFrame({
            "Total_bin_cyclo_count_Bin1_le": [10, 20],
            "cyclo_count": [100, 200],
            "Total_bin_noncyclo_count_Bin1_le": [5, 15],
            "noncyclo_count": [50, 150],
            "Total_bin_noncyclo_count_Bin2_g": [3, 7]
        })
        result = NMD_rare_steady_state_transcript(data)
        self.assertIn("test_statistic", result.columns)

class TestZScore(unittest.TestCase):
    def test_calculate_z_score(self):
        data = pd.DataFrame({
            "Isoform": ["iso1", "iso1", "iso2"],
            "test_statistic": [1.0, 2.0, 3.0]
        })
        z_scored_df = calculate_z_score(data, group_col="Isoform", stat_col="test_statistic")
        self.assertIn("z_score", z_scored_df.columns)

class TestCalculations(unittest.TestCase):
    def test_apply_hypothesis_test(self):
        data = pd.DataFrame({
            "Isoform": ["iso1", "iso1"],
            "Cyclo_TPM": [10, 20],
            "Noncyclo_TPM": [5, 15],
            "cyclo_count": [100, 200],
            "noncyclo_count": [50, 150]
        })
        result = apply_hypothesis_test(data, "Isoform", NMD_test_statistic)
        self.assertIn("test_statistic", result.columns)

    def test_process_hypothesis_test(self):
        data = pd.DataFrame({
            "Isoform": ["iso1", "iso2"],
            "Sample": ["S1", "S2"],
            "cyclo_count": [10, 20],
            "noncyclo_count": [5, 15],
            "Cyclo_TPM": [50, 100],
            "Noncyclo_TPM": [25, 50],
            "total_cyclo": [1000, 1000],
            "total_noncyclo": [500, 500]
        })
        result = process_hypothesis_test(
            filtered_data=data,
            group_col="Isoform",
            test_statistic_func=NMD_test_statistic
        )
        self.assertIn("test_statistic", result.columns)
        self.assertIn("z_score", result.columns)

class TestExpressionMatrix(unittest.TestCase):
    def test_parse_read_stats(self):
        with open("test_read_stats.txt", "w") as f:
            f.write("read1 iso1\nread2 iso1\nread3 iso2\n")
        parsed_stats = parse_read_stats("test_read_stats.txt")
        self.assertIn("iso1", parsed_stats)
        self.assertEqual(parsed_stats["iso1"].get("read1", 0), 1)

    def test_create_expression_matrix(self):
        with open("test_stats.txt", "w") as f:
            f.write("read1 iso1\nread2 iso1\nread3 iso2\n")
        matrix = create_expression_matrix("test_stats.txt")
        self.assertIsInstance(matrix, pd.DataFrame)

    def test_create_long_format(self):
        expression_matrix = pd.DataFrame({
            "S1": [10, 20],
            "S2": [30, 40]
        }, index=["iso1", "iso2"])
        sample_info = pd.DataFrame({
            "sample": ["S1", "S2"],
            "cyclo": ["cyclo", "noncyclo"],
            "haplotype": ["HP1", "HP2"]
        })
        long_format = create_long_format(expression_matrix, sample_info)
        self.assertIn("Cyclo_TPM", long_format.columns)

if __name__ == "__main__":
    unittest.main()
