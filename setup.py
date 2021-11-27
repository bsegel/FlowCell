import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="kinetics",
    version="1.0",
    author="Becca Segel",
    author_email="becca.segel@pitt.edu",
    description="Electrochemical toolkit to analyze kinetics of three electrode cell ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bsegel/FlowCell",
    packages=setuptools.find_packages()
)
