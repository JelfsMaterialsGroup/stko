import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="stko",
    version="0.0.3",
    author="Steven Bennett, Andrew Tarzia",
    author_email="s.bennett18@imperial.ac.uk",
    description="Contains molecular optimizers and property calculators for use with stk.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JelfsMaterialsGroup/stko",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
