from setuptools import setup, find_packages

setup(
    name="IsoRanker",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "pandas>=2.2.3",
        "numpy>=2.2.2",
        "matplotlib>=3.10.0",
        "statsmodels>=0.14.4",
        "seaborn>=0.13.2",
        "scipy>=1.15.1",
    ],
    entry_points={
        "console_scripts": [
            "run_analysis=IsoRanker.run_analysis:main",
        ],
    },
    author="Hank Cheng",
    description="Package for calculating test statistics, z-scores, and ranks for isoforms.",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
