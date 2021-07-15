import setuptools

setuptools.setup(
    name="stko",
    version="0.0.35",
    author="Steven Bennett, Andrew Tarzia",
    author_email="s.bennett18@imperial.ac.uk",
    description=(
        "Contains molecular optimizers and property calculators for "
        "use with stk."
    ),
    url="https://github.com/JelfsMaterialsGroup/stko",
    packages=setuptools.find_packages(),
    install_requires=(
        'scipy',
        'numpy',
        'networkx',
    ),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
