#!/usr/bin/env python

import versioneer
from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='multi_locus_analysis',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Utilities for analyzing particle trajectories',
      long_description=readme(),
      author='Bruno Beltran',
      author_email='brunobeltran0@gmail.com',
      packages=['multi_locus_analysis'],
      package_data={'multi_locus_analysis': 'vvcf_table.csv'},
      license='MIT',
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Topic :: Utilities'
      ],
      keywords='microscopy diffusion scientific',
      url='https://github.com/brunobeltran/multi_locus_analysis',
      install_requires=['numpy', 'scipy', 'pandas'],
)
