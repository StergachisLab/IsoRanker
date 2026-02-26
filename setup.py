from setuptools import setup, find_packages

setup(
    name="IsoRanker",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "matplotlib",
        "statsmodels",
        "seaborn",
        "scipy",
        "pyreadr",
        "pysam",
        "pyhpo",
        "scikit-learn"
    ],
    entry_points={
        "console_scripts": [
            "isoranker_pb_run_analysis=IsoRanker.isoranker_pb_run_analysis:main",
            "isoranker_pb_run_analysis_full=IsoRanker.isoranker_pb_run_analysis_full:main",
            "isoranker_isoquant_run_analysis=IsoRanker.isoranker_isoquant_run_analysis:main",
	        "isoranker_pb_run_analysis_minor_iso_GOE_full=IsoRanker.isoranker_pb_run_analysis_minor_iso_GOE_full:main",
            "isoranker_pb_run_analysis_full_with_isoform_allelic=IsoRanker.isoranker_pb_run_analysis_full_with_isoform_allelic:main",
            "isoranker_pb_run_analysis_full_with_isoform_allelic_omim=IsoRanker.isoranker_pb_run_analysis_full_with_isoform_allelic_omim:main",
            "isoranker_pb_run_analysis_full_with_isoform_allelic_no_novel=IsoRanker.isoranker_pb_run_analysis_full_with_isoform_allelic_no_novel:main",
            "isoranker_pb_run_analysis_full_with_isoform_allelic_protein_coding=IsoRanker.isoranker_pb_run_analysis_full_with_isoform_allelic_protein_coding:main",
            "isoranker_pb_run_analysis_full_with_isoform_allelic_hap_NMD=IsoRanker.isoranker_pb_run_analysis_full_with_isoform_allelic_hap_NMD:main"
        ],
    },
    author="Hank Cheng",
    description="Package for calculating test statistics, z-scores, and ranks for isoforms.",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
