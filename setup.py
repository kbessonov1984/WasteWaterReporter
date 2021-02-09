#!/usr/bin/env python3
import os
from distutils.core import setup
from setuptools import find_packages

author = 'Kyrylo Bessonov'


setup(
    name='wwreport',
    include_package_data=True,
    version='0.0.1',
    python_requires='>=3.7.0',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    packages=find_packages(exclude=['tests', 'databases']),
    url='https://github.com/phac-nml',
    license='GPLv3',
    author='Kyrylo Bessonov',
    author_email='kyrylo.bessonov@canada.ca',
    description=(
        'Waste Water Report Generator on Variables of Concern'),
    keywords='SARS-COV-2',
    package_dir={'wwreport': 'wwreport'},
    package_data={'mob_suite': ['config.json']},

    install_requires=[
        'pandas>=0.22.0,<=1.0.5',
        'pysam',
        'numpy<1.19'
    ],

    entry_points={
        'console_scripts': [
            'wwreport=wwreport.wwreport:main'
        ],
    },
)
