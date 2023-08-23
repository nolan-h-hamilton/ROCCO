from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as readme_file:
    long_description = readme_file.read()

setup(
    name='rocco',
    version='0.3.2',
    author='Nolan Holt Hamilton',
    author_email='nolan.hamilton@unc.edu',
    description='Robust ATAC-seq Peak Calling for Many Samples via Convex Optimization',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/nolan-h-hamilton/rocco',
    packages=find_packages(),
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    keywords='peak-caller, atac-seq, consensus-peaks',
    install_requires=[
        'numpy',
        'cvxpy',
        'pandas',
        'pararead',
        'logmuse',
        'scipy',
        'pysam',
        'pybedtools',
        'ortools',
    ],
    entry_points={
        'console_scripts': [
            'rocco = rocco.rocco:main'
        ]
    },
)