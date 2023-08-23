from setuptools import setup, find_packages

setup(
    name='rocco',
    version='0.2.0',
    packages=find_packages(),
    install_requires=['numpy==1.23.2','cvxpy==1.3.1',
                      'pandas==1.5.3', 'pararead==0.7.0',
                      'logmuse==0.2.7', 'scipy==1.10.1',
                      'pysam==0.20.0', 'pybedtools==0.9.0',
                      'ortools==9.3.10497'],
    entry_points={
        'console_scripts': [
            'rocco = rocco.rocco:main'
        ]
    }
)
