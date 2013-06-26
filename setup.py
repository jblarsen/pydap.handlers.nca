from setuptools import setup, find_packages
import sys, os


version = '0.2'

install_requires = [
    # List your project dependencies here.
    # For more details, see:
    # http://packages.python.org/distribute/setuptools.html#declaring-dependencies
    'Pydap>=3.2',
    'ConfigObj',
    'pupynere',
    'Numpy',
]


setup(name='pydap.handlers.nca',
    version=version,
    description="Handler for aggregating NetCDF data",
    long_description="",
    classifiers=[
      # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    ],
    keywords='netcdf opendap dods dap science meteorology oceanography',
    author='Roberto De Almeida',
    author_email='roberto@dealmeida.net',
    url='https://github.com/robertodealmeida/pydap.handlers.nca',
    license='MIT',
    packages=find_packages('src'),
    package_dir = {'': 'src'},
    namespace_packages = ['pydap', 'pydap.handlers'],
    include_package_data=True,
    zip_safe=False,
    install_requires=install_requires,
    entry_points="""
        [pydap.handler]
        nca = pydap.handlers.nca:NetCDFAggregatorHandler
    """,
)
