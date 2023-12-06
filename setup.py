import os
from setuptools import setup

setup(
    name = "insilicosv",
    version = "0.0.5",
    description = ("Structural variant simulation"),
    url = "https://github.com/PopicLab/insilicoSV",

    packages=['insilicosv'],
    #long_description=read('README.md'),
    long_description='Simulator of genomic structural variants that supports complex and custom variant structures',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
)
