
from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name = 'biosurfer',
    version = '0.1',
    author = 'Sheynkman Lab',
    author_email = 'sheynkman.lab@gmail.com',
    description = 'Biosurfer allows easy navigation between transcriptome and proteome',
    url = 'https://github.com/sheynkman-lab/biosurfer',
    license='MIT',
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    keywords = ['biology', 'proteomics', 'alternate splicing'],
    zip_safe = False,
    include_package_data = True,
    packages =
        ['biosurfer'],
    python_requires = '>=3.10',
    install_requires =[
        'attrs',
        'biopython < 1.87',
        'brokenaxes',
        'intervaltree >= 3.2.1',
        'matplotlib',
        'more-itertools',
        'numpy',
        'sqlalchemy >= 1.4',
        'tqdm',
        'Click',
        'seaborn',
        'ipython',
        'pandas',
        'xlsxwriter',
        'pysam',
        ],
    entry_points = {
        'console_scripts':[
            'biosurfer = biosurfer.biosurfer:cli',
        ],
    },
)
