from setuptools import setup

setup(
    name='DrillPyze',
    version='0.0.1',
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