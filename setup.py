from setuptools import setup

#
# to use this file, run:
# To check the setup file
#       python setup.py check
#
# To create the distributable
#       python setup.py sdist
#

setup(
    name='DrillPyze',
    version='0.0.2a',
    description='An accelerated drillstring dynamics simulator and estimator',
    url='https://ucalgary.ca/labs/geothermal-energy/',
    author='Roman J. Shor',
    author_email='roman.shor@ucalgary.ca',
    license='MIT',
    packages=['DrillPyze'],
    install_requires=['well_profile',
                     'numpy',
                     'pandas',
                     'pathos',
                     'scipy',
                     'numba'],
    classifiers=[
    'Development Status :: 1 - Planning',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.8',
    ]
)