SkyWalker - Python observation planner tool
-------------------------------------------

``Version: 1.0.0``

By Herpich F. R.  

This tool can be used to plan the nights for virtually any observatory on the Planet. The user can make maps for individual objects or lists containing several of them. It is also possible to define time blocks for every object (individual or in a list). The angular distance to the Moon will always be shown at the given initial time for each object (if none is given, the default is 0 LT).

Usage
+++++

- to get the full set of options available with the full description:

``python skywalker.py --help``

Requisites
++++++++++

``python 3``

``numpy``   
``matplotlib``   
``astropy``   
``datetime``   
``warnings``

This code uses a modified version of the Astroplan code (https://astroplan.readthedocs.io/en/latest/), dubbed here myastroplan. All the modified code is included in this repository. If you use this code in your research, please cite accordingly (see https://github.com/astropy/astroplan for the full reference provided by the authors).

Examples
++++++++

* Showing the track for NGC104 for Cerro Tololo and its distance to the Moon at 0:30 LT

``python skywalker.py --object NGC104 --observatory 'Cerro Tololo' -ns 2019-08-23 --at 0:30:00 --savechart test01``

.. image:: figs/test01_2019-08-23_plan.png

* Showing the skychart for the same track

``python skywalker.py --object NGC104 --observatory 'Cerro Tololo' -ns 2019-08-23 --at 0:30:00 --skychart --savechart test02``

.. image:: figs/test02_2019-08-23_plan.png
   
* Adding an observing block starting at 0:30 LT for NGC104 at Cerro Tololo

``python skywalker.py --object NGC104 --observatory 'Cerro Tololo' -ns 2019-08-23 --at 0:30:00 --blocktime 3851 --skychart --savechart test03``

.. image:: figs/test03_2019-08-23_plan.png

* Showing all tracks of a list of objects for Cerro Tololo

``python skywalker.py -f example_file.csv --observatory 'Cerro Tololo' -ns 2019-08-23 --skychart --savechart test04``

.. image:: figs/test04_2019-08-23_plan.png

* Showing all tracks of a list of objects for a given observatory provided by the sitefile

``python skywalker.py -f example_file.csv --sitefile sitefilename_example.csv -ns 2019-08-23 --skychart --savechart test05``

.. image:: figs/test05_2019-08-23_plan.png

* Including an altitude/airmass limit to the observations

``python skywalker.py -f example_file.csv --sitefile sitefilename_example.csv -ns 2019-08-23 --skychart --minalt 25 --savechart test06``

.. image:: figs/test06_2019-08-23_plan.png

Last modifications
++++++++

* 2020-03-12: Upgrading for Python 3 - Files modificated: skywalker.py and myastroplan/sky.py
