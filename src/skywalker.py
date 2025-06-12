#!/bin/python3

import os
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body
from astropy.time import Time
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import logging
from astroplan import Observer
from astroplan.plots import plot_sky
from timezonefinder import TimezoneFinder
import pytz


def parse_args():
    parser = argparse.ArgumentParser(
        description="Skywalker: A tool to visualize the sky.")
    # site
    parser.add_argument("--lat", type=float, required=False, default=-30.1652,
                        help="Latitude of the observer in degrees. Default is CTIO.")
    parser.add_argument("--lon", type=float, required=False, default=-70.8109,
                        help="Longitude of the observer in degrees. Default is CTIO.")
    parser.add_argument("--elev", type=float, required=False, default=2200,
                        help="Altitude of the observer in meters. Default is CTIO.")
    parser.add_argument("--site", type=str, required=False,
                        help="Site name (e.g., 'CTIO').")
    parser.add_argument("--sitefile", "-sf", type=str, required=False,
                        help="File containing the coordinates of the site. \
                        Can be used as an alternative to lat, lon, and alt or site.")
    parser.add_argument("--minalt", type=float, required=False, default=10,
                        help="Minimum altitude the telescope can safely go \
                        in degrees. Default is 10.")
    parser.add_argument("--sites", action="store_true",
                        help="List all available sites in the database.")

    # night
    parser.add_argument("--nightstarts", "-ns", type=str, required=False,
                        default=Time.now().strftime('%Y-%m-%d'),
                        help="Night starts (YYYY-MM-DD) in which to plot the \
                        Moon and other astronomical parameters.")
    parser.add_argument("--time", type=str, required=False, default="23:59:59",
                        help="Time (HH:MM:SS) in which to plot the Moon and other astronomical parameters.")

    # object
    parser.add_argument("--object", "-o", type=str, required=False,
                        help="Object to plot (e.g., 'M31'). Default is None.")
    parser.add_argument("--ra", required=False,
                        help="Right Ascension of the object in degrees or \
                        hourangle. If using hourangle, set --raunit to 'hour'.")
    parser.add_argument("--dec", required=False, type=str,
                        help="Declination of the object in degrees. \
                        If using hexagesimal format, use --dec='DD:MM:SS'.")
    parser.add_argument("--raunit", type=str, default='deg',
                        help="Units of the RA parameter. Options: [hour, deg]. \
                        Default is deg.")
    parser.add_argument("--obj_nme", type=str, required=False, default="Obj",
                        help="Name of the object to plot. \
                        Default is Obj.")
    parser.add_argument("--file", "-f", type=str, required=False,
                        help="File containing the coordinates of multiple objects.")
    parser.add_argument("--pid", type=str, required=False,
                        help="If file is provided and the column PID is present, \
                        only objects with the given PID will be plotted.")
    parser.add_argument("--blockinit", type=str, required=False,
                        help="Initial time of the observation block in HH:MM:SS. \
                        In case blockinit is provided, it replaces the parameter time.")
    parser.add_argument("--blocktime", type=float, required=False,
                        help="Size of the observation block in seconds.")

    # plot
    parser.add_argument("--skychart", "-sc", action="store_true",
                        help="Plot the sky chart for the target(s) starting at \
                        the time given by --time.")
    parser.add_argument("--savefig", action="store_true",
                        help="Save the figure.")
    parser.add_argument("--figname", "-fn", type=str, required=False,
                        help="Name of the figure to save. \
                        Default is skywalker_<nightstarts>.png.")

    # logging
    parser.add_argument("--logfile", type=str, default='skywalker.log',
                        help="Log file to save the log messages.")
    parser.add_argument("--loglevel", type=str, default='INFO',
                        help="Log level. Options: [DEBUG, INFO, WARNING, ERROR, CRITICAL]. \
                        Default is INFO.")

    if '-h' in os.sys.argv or '--help' in os.sys.argv:
        parser.print_help()
        import sys
        sys.exit(0)
    try:
        args = parser.parse_args()
    except:
        import sys
        e = sys.exc_info()
        parser.error(f"Argument error: {e}. \
            If using hexagesimal DEC, use --dec='DD:MM:SS' to pass the argument.")
    args = parser.parse_args()
    if args.raunit not in ['deg', 'hour']:
        parser.error("Invalid value for --raunit. Options are: [deg, hour].")
    if args.raunit == 'hour' and args.ra is None:
        parser.error(
            "If --raunit is 'hour', you must provide a value for --ra.")
    if args.raunit == 'deg' and args.ra is not None:
        try:
            float(args.ra)
        except ValueError:
            parser.error(
                "If --raunit is 'deg', --ra must be a valid float value.")
    return args


def logger(logfile=None, loglevel=logging.INFO):
    logger = logging.getLogger(__name__)

    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] @%(module)s.%(funcName)s() %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.setLevel(loglevel)

    if logfile is not None:
        fh = logging.FileHandler(logfile)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


class Skywalker:
    """
    Skywalker class to visualize the sky at a given location and time.
    It can plot the Moon, Sun, and a target object, and save the figure.
    """

    def __init__(self, args):
        """Initialize the Skywalker class with the given arguments."""

        self.location = None
        self.lat = args.lat
        self.lon = args.lon
        self.elev = args.elev
        self.site = args.site
        self.sitefile = args.sitefile
        self.sitename = None
        self.minalt = args.minalt

        self.nightstarts = args.nightstarts
        self.time = args.time if args.time else "23:59:59"
        self.inithour = None
        self.delta_midnight = None
        self.frame_time_overnight = None

        self.moon = None
        self.moon_brightness = None
        self.moonaltaz_time_overnight = None
        self.sunaltaz_time_overnight = None

        self.object = args.object
        self.ra = args.ra
        self.dec = args.dec
        self.raunit = args.raunit
        self.obj_nme = args.obj_nme
        self.file = args.file
        self.load_file = False
        self.pid = args.pid
        self.blockinit = args.blockinit
        self.blocktime = args.blocktime
        self.target = None
        self.target_list = pd.DataFrame()

        self.make_skychart = args.skychart
        self.savefig = args.savefig
        self.figname = args.figname

        self.observer = None
        self.utcoffset = -3 * u.hour

        self.logfile = args.logfile
        self.loglevel = args.loglevel
        self.logger = logger(self.logfile, self.loglevel)

    def set_location(self):
        """Set the location of the observer based on the provided parameters."""
        if self.sitefile is not None:
            if not os.path.isfile(self.sitefile):
                raise ValueError(f"File {self.sitefile} not found.")
            else:
                df = pd.read_csv(self.sitefile)
            self.location = EarthLocation(lat=df['LAT'][0],
                                          lon=df['LON'][0],
                                          height=df['ELEV'][0] * u.m)
            self.sitename = df['NAME'][0]
        elif self.site is not None:
            try:
                self.location = EarthLocation.of_site(self.site)
            except ValueError:
                raise ValueError(
                    f"Site '{self.site}' not found in the database.")
            self.sitename = self.site
        else:
            self.location = EarthLocation(lat=self.lat * u.deg,
                                          lon=self.lon * u.deg,
                                          height=self.elev * u.m)
            self.sitename = f"{self.lat} {self.lon}"
        self.logger.info(f"Location set to {self.location}")

    def set_observer(self):
        """Set the observer based on the location."""
        self.observer = Observer(self.location)

    def set_time(self):
        """Set the time of the observation."""
        if self.time is None:
            self.inithour = "23:59:59"
        else:
            if len(self.time.split(':')) != 3:
                try:
                    _time = float(self.time)
                    _time_is_str = False
                except ValueError:
                    _time_is_str = True
                if _time_is_str:
                    while len(self.time.split(':')) < 3:
                        self.time += ":00"
                else:
                    _time = f"{int(self.time)}"
                    _time += f":{int((float(self.time) - int(self.time)) * 60)}"
                    _time += ":00"
                    self.time = _time
            self.inithour = Time(self.nightstarts + "T" + self.time,
                                 format='isot').strftime('%H:%M:%S')
        self.obs_time = Time(self.nightstarts + "T" +
                             self.inithour, scale='utc', format='isot') - self.utcoffset
        tf = TimezoneFinder()
        timezone_str = tf.timezone_at(
            lng=self.location.lon.value, lat=self.location.lat.value)
        _timezone = pytz.timezone(timezone_str)
        self.utcoffset = (_timezone.utcoffset(
            self.obs_time.datetime).seconds / 3600 - 24) * u.hour

    def set_target(self, ra=None, dec=None):
        """Set the target object based on the provided parameters.
        If an object name is provided, it will try to resolve it to coordinates.
        If RA and DEC are provided, it will create a SkyCoord object.
        If a file is provided, it will load the coordinates from the file.

        Parameters:
        -----------
        ra : str or float, optional
            Right Ascension of the object in degrees or hourangle.
            If using hourangle, set --raunit to 'hour'.
        dec : str or float, optional
            Declination of the object in degrees.
        """
        if self.object:
            try:
                self.target = SkyCoord.from_name(self.object)
            except ValueError:
                raise ValueError(
                    f"Object '{self.object}' not found in the database.")
        elif (ra is not None) and (dec is not None):
            if self.raunit == 'hour':
                if (":" in ra) and (len(ra.split(':')) < 3):
                    _lenra = len(ra.split(':'))
                    while _lenra < 3:
                        ra += ":00"
                        _lenra = len(ra.split(':'))
            else:
                try:
                    ra = float(ra)
                except ValueError:
                    raise ValueError(
                        f"RA '{ra}' is not a valid float or hour format.")
                _ra = f"{int(ra)}"
                _ra += f":{int((ra - int(ra)) * 60)}"
                _ra += ":00"
                ra = _ra
            self.target = SkyCoord(ra=ra,
                                   dec=dec,
                                   unit=(self.raunit, 'deg'))
        elif self.file:
            if not os.path.isfile(self.file):
                raise ValueError(f"File {self.file} not found.")
            else:
                self.logger.info(f"File {self.file} found.")
                self.load_file = True
        else:
            raise ValueError("No target specified.")

    def set_night_frames(self):
        """Set the frames for the night observation."""
        _night_ends = (self.obs_time + .5 * u.day).strftime('%Y-%m-%d')
        _midnight = Time(f"{_night_ends}T00:00:00",
                         format='isot') - self.utcoffset
        self.delta_midnight = np.linspace(-12, 12, 500) * u.hour
        _times_time_overnight = _midnight + self.delta_midnight
        self.frame_time_overnight = AltAz(obstime=_times_time_overnight,
                                          location=self.location)
        _sun = get_body('sun', self.obs_time, self.location)
        self.moon = get_body('moon', self.obs_time, location=self.location)
        self.sunaltaz_time_overnight = get_body(
            'sun', _times_time_overnight).transform_to(self.frame_time_overnight)
        self.moonaltaz_time_overnight = get_body(
            'moon', _times_time_overnight).transform_to(self.frame_time_overnight)
        elongation = _sun.separation(self.moon)
        moon_phase = np.arctan2(_sun.distance * np.sin(elongation),
                                self.moon.distance - _sun.distance * np.cos(elongation))
        self.moon_brightness = (1. + np.cos(moon_phase)) / 2.

    def check_blockinit_format(self, blockinit=np.array([0])):
        """Check the format of the blockinit parameter and convert it to HH:MM:SS format."""
        _blockinit = blockinit.copy()
        for i, binit in enumerate(blockinit):
            try:
                _blinit = float(binit)
                _blinit_is_str = False
            except ValueError:
                _blinit_is_str = True

            if _blinit_is_str:
                while len(binit.split(':')) < 3:
                    binit += ":00"
            else:
                _blinit = f"{int(binit)}"
                _blinit += f":{int((float(binit) - int(binit)) * 60)}"
                _blinit += ":00"
                binit = _blinit
            _blockinit[i] = binit

            return _blockinit

    def set_target_list(self):
        """Set the target list based on the provided parameters."""
        if self.load_file:
            df = pd.read_csv(self.file)
            if self.pid:
                if 'PID' not in df.columns:
                    self.logger.warning(
                        f"PID column not found in {self.file}. \
                        Using all targets.")
                    self.pid = None
                else:
                    self.logger.info(
                        f"PID column found in {self.file}. \
                        Using only targets with PID {self.pid}.")
                    df = df[df['PID'] == self.pid]
            if 'RA' not in df.columns or 'DEC' not in df.columns:
                raise ValueError(
                    f"RA and DEC columns not found in {self.file}.")
            if 'NAME' not in df.columns:
                self.logger.warning(
                    f"NAME column not found in {self.file}. \
                    Using Obj incremented by 1 as names.")
                df['NAME'] = ["Obj" + str(i) for i in range(1, len(df) + 1)]
            if 'BLINIT' not in df.columns:
                self.logger.warning(
                    f"BLINIT column not found in {self.file}. \
                    Using time parameter as BLINIT.")
                df['BLINIT'] = [self.time] * len(df)
            else:
                df['BLINIT'] = self.check_blockinit_format(df['BLINIT'])
                df['BLINIT'] = df['BLINIT'].apply(
                    lambda x: Time(self.nightstarts + "T" + x, format='isot').strftime('%H:%M:%S'))
            if 'BLOCKTIME' not in df.columns:
                self.logger.warning(
                    f"BLOCKTIME column not found in {self.file}. \
                    Using 0 as BLOCKTIME.")
                df['BLOCKTIME'] = [0] * len(df)
            _coords = SkyCoord(ra=df['RA'],
                               dec=df['DEC'],
                               unit=(self.raunit, 'deg'))
            df['RA'] = _coords.ra.value
            df['DEC'] = _coords.dec.value
            self.target_list = df
        elif self.object:
            if self.blockinit:
                self.blockinit = self.check_blockinit_format(
                    np.array([self.blockinit]))
                blinit = Time(self.nightstarts + "T" + self.blockinit,
                              format='isot').strftime('%H:%M:%S')
            else:
                blinit = Time(self.nightstarts + "T" + self.time,
                              format='isot').strftime('%H:%M:%S')
            if self.blocktime:
                blocktime = self.blocktime * u.s
            else:
                blocktime = 0 * u.s
            self.obj_nme = self.object
            self.target_list = pd.DataFrame(
                {'NAME': [self.obj_nme],
                 'RA': [self.target.ra.value],
                 'DEC': [self.target.dec.value],
                 'BLINIT': [blinit],
                 'BLOCKTIME': [blocktime.value]})
        elif (self.ra is not None) and (self.dec is not None):
            if self.blockinit:
                self.blockinit = self.check_blockinit_format(
                    np.array([self.blockinit]))
                blinit = Time(self.nightstarts + "T" + self.blockinit,
                              format='isot').strftime('%H:%M:%S')
            else:
                blinit = Time(self.nightstarts + "T" + self.time,
                              format='isot').strftime('%H:%M:%S')
            if self.blocktime:
                blocktime = self.blocktime * u.s
            else:
                blocktime = 0 * u.s
            if self.obj_nme:
                obj_nme = self.obj_nme
            else:
                obj_nme = "Obj1"
            self.target_list = pd.DataFrame(
                {'NAME': [obj_nme],
                 'RA': [self.target.ra.value],
                 'DEC': [self.target.dec.value],
                 'BLINIT': [blinit],
                 'BLOCKTIME': [blocktime.value]})
        else:
            raise ValueError(
                "No target specified. Please provide a target name or coordinates.")

    def set_skychart(self,
                     observer: Observer,
                     obj_coords: SkyCoord,
                     observe_time: Time,
                     ax: plt.Axes,
                     obj_style: dict = {'color': 'b'},
                     hours_value: np.ndarray = None
                     ):
        """Set the skychart for the given object coordinates and observation time.

        Parameters:
        -----------
        observer : Observer
            The observer object containing the location and time information.
        obj_coords : SkyCoord
            The coordinates of the object to plot.
        observe_time : Time
            The time of the observation.
        ax : matplotlib.axes.Axes
            The axes on which to plot the skychart.
        obj_style : dict, optional
            Style parameters for the object plot, such as color and marker.
            Default is {'color': 'b'}.
        hours_value : np.ndarray, optional
            Array of hour values corresponding to the observation time.
            If not provided, it will be calculated from the observe_time.
        """

        try:
            plot_sky(obj_coords,
                     observer,
                     observe_time,
                     ax=ax,
                     style_kwargs=obj_style,
                     hours_value=hours_value)
        except TypeError as e:
            self.logger.error(f"Error plotting skychart: {e}")
            self.logger.error(
                "The TypeError may be due to an incompatible version of astroplan.")
            raise TypeError("Please ensure you have the modified astroplan version installed. \
                            You can get the latest version from https://github.com/herpichfr/astroplan")

    def set_plot(self):
        """Set the plot for the night observation."""
        if self.make_skychart:
            fig = plt.figure(figsize=(16, 6))
            ax1 = fig.add_subplot(121)
            ax3 = fig.add_subplot(122, projection='polar')
        else:
            fig, ax1 = plt.subplots(figsize=(8, 6))

        if self.target_list.index.size > 1:
            is_list = True
        else:
            is_list = False

        for index in self.target_list.index:
            myObjdf = self.target_list.loc[index]
            obj_coords = SkyCoord(ra=myObjdf['RA'],
                                  dec=myObjdf['DEC'],
                                  unit=('deg', 'deg'))

            block_starts = float(myObjdf['BLINIT'].split(':')[0]) + \
                float(myObjdf['BLINIT'].split(':')[1]) / 60. + \
                float(myObjdf['BLINIT'].split(':')[2]) / 3600.
            if block_starts > 12.:
                block_starts -= 24

            myaltaz_overnight = obj_coords.transform_to(
                self.frame_time_overnight)

            mask = myaltaz_overnight.alt > 0 * u.deg
            mask &= self.sunaltaz_time_overnight.alt < 0 * u.deg
            if mask.sum() == 0:
                self.logger.warning(
                    f"Object {myObjdf['NAME']} is not observable at the given \
                    time and observatory.")
                continue

            init_observable = self.frame_time_overnight[mask].obstime.min(
            )
            end_observable = self.frame_time_overnight[mask].obstime.max(
            )
            observe_time = Time(np.arange(init_observable.jd,
                                          end_observable.jd, 1./24),
                                format='jd')

            hours_values = np.array([obs_time.datetime.hour +
                                     obs_time.datetime.minute / 60.
                                     for obs_time in observe_time + self.utcoffset])

            if myObjdf['BLOCKTIME'] > 0:
                blocktime = float(myObjdf['BLOCKTIME']) * u.s
                blocktime = blocktime.to(u.hour).value
                block_ends = block_starts + blocktime
            else:
                print("Block time is 0")

            if is_list:
                p = ax1.plot(self.delta_midnight.value,
                             myaltaz_overnight.alt.value,
                             label=f"{myObjdf['NAME']}",
                             zorder=11)

                if myObjdf['BLOCKTIME'] > 0:
                    ax1.fill_between(self.delta_midnight.to('hr').value,
                                     np.zeros(
                        len(self.delta_midnight.value)),
                        myaltaz_overnight.alt.value,
                        (self.delta_midnight.value >= block_starts) & (
                        self.delta_midnight.value <= block_ends),
                        color=p[0].get_color(),
                        zorder=11)

                if self.make_skychart:
                    self.set_skychart(self.observer, obj_coords,
                                      observe_time, ax3,
                                      obj_style={'color': p[0].get_color(),
                                                 'marker': '*',
                                                 'label': myObjdf['NAME']},
                                      hours_value=hours_values)
                    print("obj_coords:", obj_coords,
                          "observe_time:", observe_time)
            else:
                sc = ax1.scatter(self.delta_midnight.value,
                                 myaltaz_overnight.alt.value,
                                 c=myaltaz_overnight.az.value,
                                 label=myObjdf['NAME'],
                                 lw=0, s=8, cmap='viridis',
                                 zorder=11)
                plt.colorbar(sc, pad=0.1).set_label('Azimuth [deg]')

                if myObjdf['BLOCKTIME'] > 0:
                    ax1.fill_between(self.delta_midnight.to('hr').value,
                                     np.zeros(
                        len(self.delta_midnight.value)),
                        myaltaz_overnight.alt.value,
                        (self.delta_midnight.value >= block_starts) & (
                        self.delta_midnight.value <= block_ends),
                        color='orange', zorder=11)

                if self.make_skychart:
                    self.set_skychart(self.observer, obj_coords,
                                      observe_time, ax3,
                                      obj_style={'cmap': 'viridis_r',
                                                 'marker': '*',
                                                 'c': hours_values,
                                                 'label': myObjdf['NAME']},
                                      hours_value=hours_values)

            # add distance to the moon at the time of the observation
            moon_distance = self.moon.separation(obj_coords).value
            text_position = abs(self.delta_midnight.value - block_starts) == abs(
                self.delta_midnight.value - block_starts).min()
            if myaltaz_overnight.alt.value[text_position].size == 0:
                self.logger.warning(
                    f"Could not find altitude for {myObjdf['NAME']} at {block_starts}")
                altitude_position = 0.0
            elif myaltaz_overnight.alt.value[text_position].size > 1:
                self.logger.warning(
                    f"Found more than one altitude for {myObjdf['NAME']} at {block_starts}")
                altitude_position = myaltaz_overnight.alt.value[text_position].mean(
                )
            else:
                altitude_position = myaltaz_overnight.alt.value[text_position]

            ax1.text(block_starts - 0.3,
                     altitude_position - 3,
                     "%i" % moon_distance,
                     fontsize=10, color='pink', zorder=12)

        ax1.plot(self.delta_midnight.to('hr').value,
                 self.moonaltaz_time_overnight.alt.value,
                 color='c', ls='--',
                 label='Moon: %i%%' % (
                     self.moon_brightness.value * 100),
                 zorder=10)
        ax1.fill_between(self.delta_midnight.to('hr').value, 0, 90,
                         (self.sunaltaz_time_overnight.alt < -0 *
                          u.deg) & (self.sunaltaz_time_overnight.alt > -18 * u.deg),
                         color='indigo', zorder=0)
        ax1.fill_between(self.delta_midnight.to('hr').value, 0, 90,
                         (self.moonaltaz_time_overnight.alt < 0 *
                          u.deg) & (self.sunaltaz_time_overnight.alt < -18 * u.deg),
                         color='k', zorder=1)
        ax1.fill_between(self.delta_midnight.to('hr').value, 0, 90,
                         (self.moonaltaz_time_overnight.alt > 0 *
                          u.deg) & (self.sunaltaz_time_overnight.alt < -18 * u.deg),
                         color='midnightblue',
                         alpha=1. - self.moon_brightness.value,
                         zorder=2)

        if self.make_skychart:
            # plot the moon into skychart
            mask = self.sunaltaz_time_overnight.alt < 0 * u.deg
            mask &= self.moonaltaz_time_overnight.alt > 0 * u.deg
            if mask.sum() > 0:
                init_moon_time = self.moonaltaz_time_overnight.obstime[mask].min(
                )
                end_moon_time = self.moonaltaz_time_overnight.obstime[mask].max(
                )
                moon_time = Time(np.arange(init_moon_time.jd,
                                           end_moon_time.jd, 1/24),
                                 format='jd')
                moon_hours = np.array([obs_time.datetime.hour +
                                       obs_time.datetime.minute / 60. for obs_time in moon_time + self.utcoffset])

                self.set_skychart(self.observer, SkyCoord(ra=self.moon.ra,
                                                          dec=self.moon.dec),
                                  moon_time.isot, ax3,
                                  obj_style={'color': 'c',
                                             'marker': 'o',
                                             'label': 'Moon: %i%%' % (self.moon_brightness.value * 100)},
                                  hours_value=moon_hours)
            else:
                self.logger.warning(
                    "Moon is not observable at the given time and observatory.")

            circle = plt.Circle((0., 0.), 90, transform=ax3.transData._b,
                                color="red", alpha=0.7, zorder=0)
            ax3.add_artist(circle)
            circle = plt.Circle((0., 0.), 90 - self.minalt, transform=ax3.transData._b,
                                color="white", alpha=1., zorder=0)
            ax3.add_artist(circle)
            circle = plt.Circle((0., 0.), 90 - self.minalt, transform=ax3.transData._b,
                                color="black", alpha=1. - self.moon_brightness.value, zorder=0)
            ax3.add_artist(circle)

        minx = self.delta_midnight.value[self.sunaltaz_time_overnight.alt < -
                                         0 * u.deg].min() - 1
        maxx = self.delta_midnight.value[self.sunaltaz_time_overnight.alt < -
                                         0 * u.deg].max() + 1
        if self.minalt > 1:
            ax1.plot([minx, maxx],
                     [self.minalt, self.minalt], '--', c='r')
        ax1.set_xlim(minx, maxx)
        ax1.set_ylim(0, 90)
        ax1.set_xlabel(f'Local Time [UTC{int(self.utcoffset.value)}]')
        xt = ax1.get_xticks()
        xt[xt < 0] += 24
        ax1.set_xticklabels(['%i' % n for n in xt])
        ax1.set_ylabel('Altitude [deg]')
        titlenight = f"Night starts: {self.nightstarts} @ {self.sitename}"
        ax1.set_title(titlenight, fontsize=11)
        ax1.legend(loc='upper right', fontsize=8)

        ax2 = ax1.twinx()
        altitude = ax1.get_yticks() * u.deg
        airmass = 1. / np.cos(90 * u.deg - altitude)
        ax2.set_ylabel('Airmass')
        myticks = []
        for airval in airmass:
            if airval > 10:
                myticks.append('')
            else:
                myticks.append('%.2f' % airval)
        ax2.set_yticklabels(myticks)
        ax2.set_ylim(0, 90)

        if self.make_skychart:
            ax3.legend(loc='lower right', fontsize=8,
                       bbox_to_anchor=(1., -0.1))

        plt.grid()
        plt.tight_layout()

        if self.savefig:
            if self.figname:
                fig.savefig(self.figname, dpi=300,
                            bbox_inches='tight')
                self.logger.info(
                    f"Figure saved as {self.figname}")
                print(f"Figure saved as {self.figname}")
            else:
                fig.savefig(f"skywalker_{self.nightstarts}.png",
                            dpi=300, bbox_inches='tight')
                self.logger.info(
                    f"Figure saved as skywalker_{self.nightstarts}.png")
                print(f"Figure saved as skywalker_{self.nightstarts}.png")
        plt.show()

    def main(self):
        self.set_location()
        self.set_observer()
        self.set_time()
        self.set_target(ra=self.ra, dec=self.dec)
        self.set_night_frames()
        self.set_target_list()
        self.set_plot()


if __name__ == "__main__":
    args = parse_args()
    if args.sites:
        print("Available sites:")
        for site in EarthLocation.get_site_names():
            print(site)
    else:
        Luke = Skywalker(args)
        Luke.main()
