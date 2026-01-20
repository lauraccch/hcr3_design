from setuptools import setup, find_packages

setup(
    name="hcr3_design",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy==2.3.2",
        "pandas==2.3.2",
        "biopython==1.85"
    ],
    python_requires=">=3.13",
    author="Laura Cecchino",
    description="HCR3 probe design",
)
