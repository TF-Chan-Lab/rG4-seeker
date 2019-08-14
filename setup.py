"""
**************************************************************************
**  rG4-seeker
**  A pipeline for processing and analyzing rG4-seq data
**
**  Version 1.0 -- August 10, 2019
**
**  Copyright (C) 2019 by Eugene Yui-Ching Chow, Ting-Fung Chan, All rights reserved.
**  Contact:  eugenechow823@gmail.com, tf.chan@cuhk.edu.hk
**  Organization:  School of Life Sciences, The Chinese University of Hong Kong,
**                 Shatin, NT, Hong Kong SAR
**
**  This file is part of rG4-seeker.
**
**  rG4-seeker is free software; you can redistribute it and/or
**  modify it under the terms of the GNU General Public License
**  as published by the Free Software Foundation; either version
**  3 of the License, or (at your option) any later version.
**
**  rG4-seeker is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public
**  License along with rG4-seeker; if not, see
**  <http://www.gnu.org/licenses/>.
**************************************************************************/
"""

from setuptools import setup, find_packages
from rg4seeker import __version__, __author__, __email__

with open('requirements.txt') as f:
    requirements = [line for line in f.read().splitlines() if line]

setup(
    name='rG4-seeker',
    version=__version__,
    packages=find_packages(),
    author=__author__,
    author_email=__email__,
    keywords='rG4, rG4-seq',
    description='rG4-seeker, A pipeline for processing and analyzing rG4-seq data',
    url='https://github.com/TF-Chan-Lab/rG4-seeker',
    include_package_data=True,
    zip_safe=False,

    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'rG4-seeker = rg4seeker.__main__:main',
        ]
    },
license = 'GPL3.0',
)
