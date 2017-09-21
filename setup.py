""" setup freeCoG for installation """

from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import io
import codecs
import os
import sys

import freeCoG

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name='freeCoG',
    version=sandman.__version__,
    url='https://github.com/ZacharyIG/freeCoG/tree/master',
    license='Apache Software License',
    author='Zachary GReenberg',
    tests_require=['pytest'],
    install_requires=['freesurfer>=0.5.3',
                    'fsl>=5.0',
                    'Python==0.2.7',
                    ],
    cmdclass={'test': PyTest},
    author_email='',
    description='Automated tools for Brain Imaging and Neuroscience Data',
    long_description=long_description,
    packages=['freeCoG'],
    include_package_data=True,
    platforms='any',
    test_suite='freeCoG.test.test_freeCoG',
    classifiers = [
        'Programming Language :: Python',
        'Development Status :: 2 - Beta',
        'Natural Language :: English',
        'Environment :: Unix/Linux',
        'Intended Audience :: Scientists-Clinicians',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Software Development :: Libraries :: Application Frameworks',
        'Topic :: Internet :: WWW/HTTP :: Dynamic Content',
        ],
    extras_require={
        'testing': ['pytest'],
    }
)
