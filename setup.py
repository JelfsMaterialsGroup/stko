import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="stko",
    version="0.0.1",
    author="Steven Bennett, Andrew Tarzia",
    author_email="s.bennett18@imperial.ac.uk",
    description="Contains optimizers for use with stk.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/stevenbennett96/stk_optimizers-1",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
