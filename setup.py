#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools.command.test import test as TestCommand

import io
import sys

import re
VERSIONFILE = "co2sys/__init__.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt', 'CHANGES.txt')


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['--strict', '--verbose',
                          '--tb=long', 'seawater/test']
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

source = 'http://pypi.python.org/packages/source/c/co2sys'
download_url = '%s/co2sys-%s.tar.gz' % (source, verstr)


README = open('README.txt').read()
CHANGES = open('CHANGES.txt').read()
LICENSE = open('LICENSE.txt').read()

setup(
    name='co2sys',
    version=verstr,
    url='http://pypi.python.org/pypi/co2sys/',
    license='MIT',
    author='Filipe Fernandes',
    tests_require=['pytest'],
    install_requires=['numpy'],
    extras_require={'testing': ['pytest']},
    cmdclass={'test': PyTest},
    author_email='ocefpaf@gmail.com',
    description='CO2SYS Library for Python',
    long_description=long_description,
    packages=['co2sys', 'co2sys.test'],
    include_package_data=True,
    platforms='any',
    test_suite='co2sys.test',
    classifiers=['Development Status :: 1 - Planning',
                 'Programming Language :: Python',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: MIT License',
                 'Operating System :: OS Independent',
                 'Topic :: Scientific/Engineering',
                 ],
    use_2to3=True,
    maintainer='Filipe Fernandes',
    maintainer_email='ocefpaf@gmail.com',
    download_url=download_url,
    keywords=['oceanography', 'co2sys', 'carbonate'],
    )
