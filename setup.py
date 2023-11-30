import os
from setuptools import setup

setup(
    name = "insilicosv",
    version = "0.0.4",
    description = ("Simulator of complex structural variants in the genome "),
    #license = "BSD",
    #keywords = "DNA simulation structural variant",
    url = "https://github.com/PopicLab/insilicoSV",

    packages=['insilicosv'],
    #long_description=read('README.md'),
    long_description='Simulator of complex genomic structural variants that support custom variant structures',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
