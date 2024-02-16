from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as readme_file:
    long_description = readme_file.read()

core_dependencies = [
    'cvxpy',
    'numpy',
    'scipy',
    'pandas',
    'pybedtools',
    'matplotlib',
    'pyBigWig',
    'deeptools'
]

optional_feature_dependencies = {
    'mosek': ['mosek'],
    'ortools': ['ortools'],
    'pytest': ['pytest']
}

optional_dependencies_message = (
    "Additional dependencies for optional features:\n\n"
    "- 'mosek': Commercial grade solver. Users can instantly obtain a free academic license or generous trial commericial license at https://www.mosek.com/products/academic-licenses/.\n"
    "- 'ortools': includes the first-order solver, PDLP.\n"
    "- 'pytest': allows local execution of the Tests workflow.\n"
)

long_description += "\n\n" + optional_dependencies_message

setup(
    name='rocco',
    version='0.11.1',
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
    include_package_data=True,
    keywords='peak-caller, atac-seq, consensus-peaks',
    install_requires=core_dependencies,
    extras_require=optional_feature_dependencies,
    entry_points={
        'console_scripts': [
            'rocco = rocco.rocco:main'
        ]
    },
)
