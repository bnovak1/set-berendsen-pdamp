from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as readme_file:
    readme = readme_file.read()

setup(
    name="set-berendsen-pdamp",
    version="0.1.0",
    packages=find_packages(include=["set_berendsen_pdamp", "set_berendsen_pdamp.*"]),
    description="Set Pdamp for Berendsen barostat to obtain a target relaxation time.",
    keywords="Berendsen barostat, LAMMPS, Pdamp, optimization",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Brian Novak",
    author_email="tajrkala@gmail.com",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8, <3.12",
    install_requires=["lmfit", "numpy", "pytest", "pandas", "scipy", "mpi4py", "openmpi"],
    test_suite="tests",
    project_urls={
        "Source": "https://github.com/username/set_berendsen_pdamp",
    },
)
