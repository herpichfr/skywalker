SkyWalker - Python observation planner tool
===========================================

[![Version](https://img.shields.io/github/v/release/herpichfr/skywalker)](https://img.shields.io/github/v/release/herpichfr/skywalker)
![GitHub issues](https://img.shields.io/github/issues/herpichfr/skywalker)
[![License](https://img.shields.io/badge/license-GNUv3.0-green)](LICENSE)
[![Python 3](https://img.shields.io/badge/python-3.6%2B-blue.svg)](https://www.python.org/downloads/)
![GitHub](https://img.shields.io/github/stars/herpichfr/skywalker?style=social)

By Herpich F. R.  

This tool can be used to plan the nights for virtually any observatory on the Planet. The user can make maps for individual objects or lists containing several of them. It is also possible to define time blocks for every object (individual or in a list). The angular distance to the Moon will always be shown at the given initial time for each object (if none is given, the default is 0 LT).

Usage
+++++

- to get the full set of options available with the full description:

``python skywalker.py --help``

Requisites
++++++++++

``python 3``

``os``
``pandas``
``timezonefinder``
``pytz``
``numpy``
``matplotlib``
``astropy``
``astroplan``

This code uses a modified version of the Astroplan code (https://astroplan.readthedocs.io/en/latest/). If you use this code in your research, please cite accordingly (see https://github.com/astropy/astroplan for the full reference provided by the authors).

Installation
++++++++++++

The package was only tested on Python 3.6 and above on Linux systyems.

To check for dependencies, run:

``bash install.sh --check``

This will check if all the required packages are installed. If any of them is missing, it will print a message with the name of the package and how to install it.

To install the required packages and dependencies, run (do not run the command with sudo):

``bash install.sh --install``

This will create a python virtual environment and install all the required packages in it, activating it if the install is successful. 

To uninstall the package, run:

``bash install.sh --uninstall``

Usage examples
++++++++

* Showing the track for NGC104 for Cerro Tololo and its distance to the Moon at 0:30 LT

``python skywalker.py --object NGC104 --site 'Cerro Tololo' -ns 2019-08-23 --at 0:30:00 --savefig test01``

.. image:: figs/test01_2019-08-23_plan.png

* Showing the skychart for the same track

``python skywalker.py --object NGC104 --site 'Cerro Tololo' -ns 2019-08-23 --at 0:30:00 --skychart --savefig test02``

.. image:: figs/test02_2019-08-23_plan.png
   
* Adding an observing block starting at 0:30 LT for NGC104 at Cerro Tololo

``python skywalker.py --object NGC104 --site 'Cerro Tololo' -ns 2019-08-23 --at 0:30:00 --blocktime 3851 --skychart --savefig test03``

.. image:: figs/test03_2019-08-23_plan.png

* Showing all tracks of a list of objects for Cerro Tololo

``python skywalker.py -f example_file.csv --site 'Cerro Tololo' -ns 2019-08-23 --skychart --savefig test04``

.. image:: figs/test04_2019-08-23_plan.png

* Showing all tracks of a list of objects for a given observatory provided by the sitefile

``python skywalker.py -f example_file.csv --sitefile sitefilename_example.csv -ns 2019-08-23 --skychart --savefig test05``

.. image:: figs/test05_2019-08-23_plan.png

* Including an altitude/airmass limit to the observations

``python skywalker.py -f example_file.csv --sitefile sitefilename_example.csv -ns 2019-08-23 --skychart --minalt 25 --savefig test06``

.. image:: figs/test06_2019-08-23_plan.png

Last modifications
++++++++

* 2020-03-12: Upgrading for Python 3 - Files modificated: skywalker.py and myastroplan/sky.py
* 2025-05-21: Full refactoring of skywalker.py. Removing myastroplan in favour of a locally modified fork of astroplan. The modified version can be found in (https://github.com/herpichfr/astroplan).
