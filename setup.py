from setuptools import setup, find_packages
import sys, os

version = '0.1dev'

setup(name                 = 'gosr',
      version              = version,
      description          = "Tools and packages for dealing with NGS data processing",
      long_description     = "",
      classifiers          = [], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords             = '',
      author               = 'Wolfgang Resch',
      author_email         = 'wresch@mail.nih.gov',
      url                  = '',
      license              = 'MIT',
      packages             = find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts              = ['bin/gosr'],
      include_package_data = True,
      zip_safe             = False,
      install_requires     = [
          # -*- Extra requirements: -*-
      ],
      entry_points         = """
      # -*- Entry points: -*-
      """,
      )
