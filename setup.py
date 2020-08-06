#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()


setup(
    name='scbrowse',
    use_scm_version={
        'local_scheme': 'dirty-tag',
        'write_to': 'src/scbrowse/_version.py',
        'fallback_version': '0.0.0',
    },
    license='LGPL-3.0-or-later',
    description='Interactive browser for single-cell ATAC-seq data',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
    ),
    author='Wolfgang Kopp',
    author_email='wolfgang.kopp@mdc-berlin.de',
    url='https://github.com/wkopp/python-scbrowse',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)'
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Utilities',
        'Private :: Do Not Upload',
    ],
    project_urls={
        'Documentation': 'https://python-scbrowse.readthedocs.io/',
        'Changelog': 'https://python-scbrowse.readthedocs.io/en/latest/changelog.html',
        'Issue Tracker': 'https://github.com/wkopp/python-scbrowse/issues',
    },
    keywords=[
        'single-cell ATAC-seq', 'interactive', 'browser',
    ],
    python_requires='>=3.7',
    install_requires=[
        'dash>=1.13',
        'dash-core-components>=1.10',
        'dash-html-components>=1.0',
        'plotly>=4.9',
        'scipy>=1.5',
        'numpy',
        'pandas',
        'pybedtools>=0.8',
        'flask_caching',
        'scregseg',
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    setup_requires=[
        'setuptools_scm>=3.3.1',
    ],
    entry_points={
        'console_scripts': [
            'scbrowse = scbrowse.cli:main',
        ]
    },
)
