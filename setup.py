from insilicosv import __version__
from setuptools import setup

setup(
    name = "insilicosv",
    version = __version__,
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
