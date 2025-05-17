#!/bin/python3

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
    parser.add_argument("--site", type=str, required=False, default=None,
                        help="Site name (e.g., 'CTIO').")
    parser.add_argument("--sitefile", "-sf", type=str, required=False,
                        help="File containing the coordinates of the site. \
                        Can be used as an alternative to lat, lon, and alt or site.")

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
    parser.add_argument("--ra", type=float, required=False,
                        help="Right Ascension of the object in degrees or hourangle.")
    parser.add_argument("--dec", type=float, required=False,
                        help="Declination of the object in degrees.")
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
                        help="Initial time of the observation block in HH:MM:SS.")
    parser.add_argument("--blocktime", type=float, required=False,
                        help="Size of the observation block in seconds.")

    # observatory
    parser.add_argument("--observatory", "-obs", type=str, required=False,
                        help="Observatory name (e.g., 'CTIO'). Default is None.")
    parser.add_argument("--observatoryfile", "-of", type=str, required=False,
                        help="File containing the coordinates of the observatory. \
                        Can be used as an alternative to observatory.")
    parser.add_argument("--minalt", type=float, required=False, default=30,
                        help="Minimum altitude of the telescope in degrees. Default is 30.")

    # plot
    parser.add_argument("--skychart", "-sc", action="store_true",
                        help="Plot the sky chart for the target(s) starting at \
                        the time given by --time.")
    parser.add_argument("--noshow", action="store_true",
                        help="Do not show the plot.")

    return parser.parse_args()


def logger():
    logger = logging.getLogger("skywalker")
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


class Skywalker:
    def __init__(self, args):
        self.location = None
        self.lat = args.lat
        self.lon = args.lon
        self.elev = args.elev
        self.site = args.site
        self.sitefile = args.sitefile

        self.nightstarts = args.nightstarts
        self.time = args.time if args.time else "23:59:59"
        self.inithour = None
        self.delta_midnight = None
        self.frame_tonight = None
        self.frame_time_overnight = None

        self.sun = None
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
        self.pid = args.pid
        self.blockinit = args.blockinit
        self.blocktime = args.blocktime
        self.target = None
        self.target_list = []

        self.observatory = args.observatory
        self.observatoryfile = args.observatoryfile
        self.minalt = args.minalt

        self.skychart = args.skychart
        self.noshow = args.noshow

        self.timezone = None
        self.nighttime = None
        self.observer = None
        self.utcoffset = None

    def set_location(self):
        if self.sitefile:
            self.location = EarthLocation.from_file(self.sitefile)
        elif self.site:
            self.location = EarthLocation.of_site(self.site)
        else:
            self.location = EarthLocation(lat=self.lat * u.deg,
                                          lon=self.lon * u.deg,
                                          height=self.elev * u.m)
        logger().info(f"Location set to {self.location}")

    def set_observer(self):
        self.observer = Observer(self.location)

    def set_time(self):
        self.nighttime = Time(self.nightstarts + "T" +
                              self.time, scale='utc', format='isot')
        logger().info(f"Time set to {self.time}")
        if not self.time:
            self.inithour = "23:59:59"
        else:
            self.inithour = Time(self.nightstarts + "T" + self.time,
                                 format='isot').strftime('%H:%M:%S')
        tf = TimezoneFinder()
        timezone_str = tf.timezone_at(lng=self.lon, lat=self.lat)
        self.timezone = pytz.timezone(timezone_str)
        dt = self.nighttime.datetime
        self.utcoffset = (self.timezone.utcoffset(
            dt).seconds / 3600 - 24) * u.hour
        self.obs_time = Time(self.nightstarts + "T" +
                             self.inithour, scale='utc', format='isot') - self.utcoffset

    def set_target(self, ra=None, dec=None):
        if self.object:
            self.target = SkyCoord.from_name(self.object)
        elif ra and dec:
            self.target = SkyCoord(ra=self.ra * u.deg,
                                   dec=self.dec * u.deg)
        else:
            raise ValueError("No target specified.")

    def set_night_frames(self):
        observer_time = self.nighttime + self.utcoffset
        night_ends = observer_time + 1 * u.day
        midnight = Time(night_ends.jd - 0.5, format='jd') + self.utcoffset
        self.delta_midnight = np.linspace(-12, 12, 500) * u.hour
        self.frame_tonight = AltAz(obstime=observer_time + self.delta_midnight,
                                   location=self.location)
        times_time_overnight = midnight + self.delta_midnight
        self.frame_time_overnight = AltAz(obstime=times_time_overnight,
                                          location=self.location)
        self.sun = get_body('sun', self.obs_time, self.location)
        self.moon = get_body('moon', self.obs_time, location=self.location)
        self.sunaltaz_time_overnight = get_body(
            'sun', times_time_overnight).transform_to(self.frame_time_overnight)
        self.moonaltaz_time_overnight = get_body(
            'moon', times_time_overnight).transform_to(self.frame_time_overnight)
        elongation = self.sun.separation(self.moon)
        moon_phase = np.arctan2(self.sun.distance * np.sin(elongation),
                                self.moon.distance - self.sun.distance * np.cos(elongation))
        self.moon_brightness = (1. + np.cos(moon_phase)) / 2.

    def set_target_list(self):
        if self.file:
            df = pd.read_csv(self.file)
            if self.pid:
                df = df[df['PID'] == self.pid]
            self.target_list = df
        elif self.object:
            try:
                target = SkyCoord.from_name(self.object, frame='icrs')
            except ValueError:
                raise ValueError(
                    f"Object '{self.object}' not found in the database.")
            if self.blockinit:
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
                 'RA': [target.ra.value],
                 'DEC': [target.dec.value],
                 'BLINIT': [blinit],
                 'BLOCKTIME': [blocktime.value]})
        elif self.ra and self.dec:
            self.target = SkyCoord(ra=self.ra,
                                   dec=self.dec, unit=(args.raunit, 'deg'))
            if self.blockinit:
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

    def skychart(self,
                 ax,
                 mycoords,
                 observe_time,
                 obj_style={'color': 'b'},
                 hours_value=False
                 ):
        plot_sky(mycoords,
                 self.observer,
                 observe_time,
                 ax=ax,
                 style_kwrds=obj_style,
                 hours_value=hours_value)

    def set_plot(self):
        if args.skychart:
            fig = plt.figure(figzise=(16, 6))
            ax1 = fig.add_subplot(121)
            ax3 = fig.add_subplot(122, projection='polar')
        else:
            fig, ax1 = plt.subplots(figsize=(16, 6))

        for myObj in self.target_list:
            objtime = Time(self.nightstarts + "T" + self.time,
                           format='isot') - self.utcoffset
            mycoords = SkyCoord(ra=myObj.ra, dec=myObj.dec,
                                obstime=objtime, frame='icrs')
            # convert self.time to decimal hours
            obj_dec_inihour = self.time.decimalhour
            inivalue = self.time.decimalhour
            endvalue = delta_midnight[self.sunaltaz_time_overnight.alt < -
                                      18 * u.deg].max().value
            observe_time = self.obs_time + \
                np.arange(inivalue, endvalue, 1) * u.hour

            myaltaz = myObj.transform_to(AltAz(obstime=observe_time,
                                               location=self.location))
            myaltaz_tonight = myObj.transform_to(self.frame_tonight)
            myaltaz_overnight = myObj.transform_to(
                self.frame_time_overnight)

            p = ax1.plot(self.delta_midnight.value,
                         myaltaz_overnight.alt.value, label=myObj['name'])

            if args.skychart:
                self.skychart(ax3, self.location, mycoords, self.obs_time,
                              observe_time,
                              self.sunaltaz_time_overnight.alt.value,
                              obs_style={'color': p[0].get_color(),
                                         'marker': '*',
                                         'label': myObj['name']},
                              hours_value=np.arange(inivalue, endvalue, 1))
            else:
                sc = ax1.scatter(self.delta_midnight.value,
                                 myaltaz_overnight.alt.value,
                                 c=myaltaz_overnight.az.value,
                                 label=myObj['name'],
                                 lw=0, s=8, cmap='viridis')
                if args.skychart:
                    self.skychart(ax3, self.location, mycoords,
                                  self.obs_time + self.utcoffset, observe_time,
                                  self.sunaltaz_time_overnight.alt.value,
                                  obs_style={'cmap': 'viridis',
                                             'marker': '*',
                                             'c': np.arange(inivalue, endvalue, 1),
                                             'label': myObj['name']},
                                  hours_value=np.arange(inivalue, endvalue, 1))

        ax1.plot(self.delta_midnight.value,
                 self.moonaltaz_time_overnight.alt.value,
                 color='violet', ls='--', label='Moon')
        ax1.fill_between(self.delta_midnight.to('hr').value, 0, 90,
                         (self.sunaltaz_time_overnight.alt < -0 *
                          u.deg) & (self.sunaltaz_time_overnight.alt > -18 * u.deg),
                         color='indigo', zorder=0)
        ax1.fill_between(self.delta_midnight.to('hr').value, 0, 90,
                         (self.moon_altaz_time_overnight.alt < 0 *
                          u.deg) & (self.moonaltaz_time_overnight.alt < -18 * u.deg),
                         color='k', zorder=0)
        ax1.fill_between(self.delta_midnight.to('hr').value, 0, 90,
                         (self.moonaltaz_time_overnight.alt > 0 *
                          u.deg) & (self.moonaltaz_time_overnight.alt < -18 * u.deg),
                         color='midnightblue', alpha=1.1 - moon_brightness, zorder=0)

    def main(self):
        self.set_location()
        self.set_observer()
        self.set_time()
        self.set_target()
        self.set_night_frames()
        self.set_target_list()


if __name__ == "__main__":
    args = parse_args()
    Luke = Skywalker(args)
    Luke.main()
