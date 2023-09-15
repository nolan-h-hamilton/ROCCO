from setuptools import setup, find_packages

with open("README_PyPI.md", "r", encoding="utf-8") as readme_file:
    long_description = readme_file.read()

# Core dependencies
core_dependencies = [
    'numpy',
    'cvxpy',
    'pandas',
    'scipy',
    'pysam',
    'pararead',
    'logmuse'
]

optional_feature_dependencies = {
    'mosek': ['mosek'],
    'ortools': ['ortools'],
}

all_dependencies = core_dependencies + sum(optional_feature_dependencies.values(), [])

optional_dependencies_message = (
    "Additional dependencies for optional features:\n\n"
    "- 'mosek': Commercial grade solver. Users can instantly obtain a free academic license or generous trial commericial license at https://www.mosek.com/products/academic-licenses/.\n"
    "- 'ortools': Includes the PDLP first-order solver.\n"
)

long_description += "\n\n" + optional_dependencies_message

setup(
    name='rocco',
    version='0.4.2',
    author='Nolan Holt Hamilton',
    author_email='nolan.hamilton@unc.edu',
    description='Robust ATAC-seq Peak Calling for Many Samples via Convex Optimization',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/nolan-h-hamilton/rocco',
    packages=find_packages(),
    scripts=
        ['rocco/smoothWig.pl', 'rocco/cutsToWig.pl'],
    include_package_data=True,
    package_data={'': ['rocco/cutsToWig.pl', 'rocco/smoothWig.pl']},
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
    ],
    keywords='peak-caller, atac-seq, consensus-peaks',
    install_requires=core_dependencies,
    extras_require=optional_feature_dependencies,
    entry_points={
        'console_scripts': [
            'rocco = rocco.rocco:main'
        ]
    },
)
