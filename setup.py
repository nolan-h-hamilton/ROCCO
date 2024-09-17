from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as readme_file:
    long_description = readme_file.read()

setup(
    name='rocco',
    version='1.1.0',
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
        "Development Status :: 4 - Beta",
    ],
    python_requires='>=3.8,<3.12',
    install_requires=[
        'ortools>=9.10',
        'numpy>=1.23.5',
        'scipy>=1.11',
        'pandas>=2.0',
        'pybedtools>=0.9',
        'pyBigWig>=0.3.22',
        'deeptools>=3.5',
        'pysam>=0.20.0',
        'myst-parser>=0.15.0',
    ],
    extras_require={
        'pytest': ['pytest>=7.0.0']
    },
    entry_points={
        'console_scripts': [
            'rocco = rocco.rocco:main'
        ]
    },
    include_package_data=True,
)
