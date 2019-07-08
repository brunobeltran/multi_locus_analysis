#!/usr/bin/env python

import versioneer
from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='multi_locus_analysis',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Utilities for analyzing particle trajectories',
      long_description=readme(),
      author='Bruno Beltran',
      author_email='brunobeltran0@gmail.com',
      packages=['multi_locus_analysis', 'multi_locus_analysis.examples',
                # 'multi_locus_analysis.examples.garcia',
                'multi_locus_analysis.examples.burgess'],
      package_data={
          'multi_locus_analysis':
              ['vvcf_table.csv'],
          # 'multi_locus_analysis.examples.burgess':
          #     ['examples/burgess/xyz_conf_okaycells9exp.csv'],
      },
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
      install_requires=['bruno_util', 'matplotlib', 'numpy', 'scipy', 'pandas',
                        'seaborn', 'mpmath', 'statsmodels'],
)
