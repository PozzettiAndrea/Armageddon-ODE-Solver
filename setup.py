#!/usr/bin/env python

from setuptools import setup, Extension

setup(name='Armageddon',
      version='1.0',
      description='Asteroid atmospheric entry solver',
      author='AMCG project',
      packages=['armageddon_model', 'acse4p2score', 'acse4p2testdata'],
      package_dir = {'armageddon_model': 'armageddon',
                     'acse4p2testdata': 'test-data'},
      package_data = {'armageddon_model':['resources/*.csv'],
                      'acse4p2score': ['pytest.ini'],
                      'acse4p2testdata':['*.csv']},
      include_package_data=True
     )
