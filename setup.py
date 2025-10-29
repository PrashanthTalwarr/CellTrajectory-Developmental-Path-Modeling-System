from setuptools import setup, find_packages

setup(
    name="celltrajectory",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "scanpy>=1.8.1",
        "anndata>=0.7.6",
        "scvelo>=0.2.4",
        "cellrank>=1.5.1",
        "scikit-learn>=1.0.0",
        "pot>=0.8.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "tqdm>=4.61.0",
        "h5py>=3.3.0",
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A framework for modeling cellular differentiation trajectories",
    keywords="single-cell, RNA-seq, trajectory inference, developmental biology",
    url="https://github.com/yourusername/CellTrajectory",
)