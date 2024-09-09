from setuptools import setup, find_packages

setup(
    name="Clust_GO",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        "requests",
        "goatools",
        "pandas",
        "numpy",
        "scikit-learn",
        "matplotlib",
        "bioservices"
    ],
    description="A package for protein analysis and GO term association",
    author="Karen Sobriel; Grégoird Ménard",
    author_email="crcina.impact.bioinfo@gmail.com",
    url="https://gitlab.univ-nantes.fr/E179974Z/clust_go",
)
