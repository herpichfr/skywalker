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

    return parser.parse_args()


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
    def __init__(self, args):
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
        self.observer = Observer(self.location)

    def set_time(self):
        if not self.time:
            self.inithour = "23:59:59"
        else:
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
        if self.object:
            try:
                self.target = SkyCoord.from_name(self.object)
            except ValueError:
                raise ValueError(
                    f"Object '{self.object}' not found in the database.")
        elif ra and dec:
            self.target = SkyCoord(ra=self.ra,
                                   dec=self.dec,
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

    def check_blockinit_format(self):
        if len(self.blockinit.split(':')) != 3:
            _blockinit = self.blockinit.split(':')
            if len(_blockinit) == 2:
                self.blockinit = f"{_blockinit[0]}:{_blockinit[1]}:00"
            elif len(_blockinit) == 1:
                self.blockinit = f"{_blockinit[0]}:00:00"

    def set_target_list(self):
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
            if 'BLOCKTIME' not in df.columns:
                self.logger.warning(
                    f"BLOCKTIME column not found in {self.file}. \
                    Using 0 as BLOCKTIME.")
                df['BLOCKTIME'] = [0] * len(df)
            self.target_list = df
        elif self.object:
            if self.blockinit:
                self.check_blockinit_format()
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
        elif self.ra and self.dec:
            if self.blockinit:
                self.check_blockinit_format()
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
                     observer,
                     obj_coords,
                     observe_time,
                     ax,
                     obj_style={'color': 'b'},
                     hours_value=None
                     ):

        plot_sky(obj_coords,
                 observer,
                 observe_time,
                 ax=ax,
                 style_kwargs=obj_style,
                 hours_value=hours_value)

    def set_plot(self):
        if self.make_skychart:
            fig = plt.figure(figsize=(16, 6))
            ax1 = fig.add_subplot(121)
            ax3 = fig.add_subplot(122, projection='polar')
        else:
            fig, ax1 = plt.subplots(figsize=(8, 6))
        plt.grid()

        if self.target_list.index.size > 1:
            is_list = True
        else:
            is_list = False

        for index in self.target_list.index:
            myObjdf = self.target_list.loc[index]
            obj_coords = SkyCoord(ra=myObjdf['RA'],
                                  dec=myObjdf['DEC'],
                                  unit=(self.raunit, 'deg'))

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
                                          end_observable.jd, 1./24), format='jd')

            hours_values = np.array([obs_time.datetime.hour + self.utcoffset.value +
                                     obs_time.datetime.minute / 60. for obs_time in observe_time])

            if myObjdf['BLOCKTIME'] > 0:
                blocktime = float(myObjdf['BLOCKTIME']) * u.s
                blocktime = blocktime.to(u.hour).value
                block_ends = block_starts + blocktime
            else:
                print("Block time is 0")

            if is_list:
                p = ax1.plot(self.delta_midnight.value,
                             myaltaz_overnight.alt.value,
                             label=myObjdf['NAME'],
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
                                      obj_style={'cmap': 'viridis',
                                                 'marker': '*',
                                                 'c': hours_values,
                                                 'label': myObjdf['NAME']},
                                      hours_value=hours_values)

            # add distance to the moon at the time of the observation
            moon_distance = self.moon.separation(obj_coords).value
            text_position = abs(self.delta_midnight.value - block_starts) == \
                abs(self.delta_midnight.value - block_starts).min()
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
                 label='Moon: %i' % (self.moon_brightness.value * 100),
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
                         alpha=1.1 - self.moon_brightness.value,
                         zorder=2)

        if self.make_skychart:
            # plot the moon into skychart
            mask = self.sunaltaz_time_overnight.alt < 0 * u.deg
            mask &= self.moonaltaz_time_overnight.alt > 0 * u.deg
            init_moon_time = self.moonaltaz_time_overnight.obstime[mask].min()
            end_moon_time = self.moonaltaz_time_overnight.obstime[mask].max()

            moon_time = Time(np.arange(init_moon_time.jd,
                                       end_moon_time.jd, 1/24), format='jd')
            moon_hours = np.array([obs_time.datetime.hour + self.utcoffset.value +
                                   obs_time.datetime.minute / 60. for obs_time in moon_time])

            self.set_skychart(self.observer, SkyCoord(ra=self.moon.ra,
                                                      dec=self.moon.dec),
                              moon_time.isot, ax3,
                              obj_style={'color': 'c',
                                         'marker': 'o',
                                         'label': 'Moon: %i' % (self.moon_brightness.value * 100)},
                              hours_value=moon_hours)

            circle = plt.Circle((0., 0.), 90, transform=ax3.transData._b,
                                color="red", alpha=0.7, zorder=0)
            ax3.add_artist(circle)
            circle = plt.Circle((0., 0.), 90 - self.minalt, transform=ax3.transData._b,
                                color="white", alpha=1., zorder=0)
            ax3.add_artist(circle)
            circle = plt.Circle((0., 0.), 90 - self.minalt, transform=ax3.transData._b,
                                color="black", alpha=1.1 - self.moon_brightness.value, zorder=0)
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

        ax3.legend(loc='lower right', fontsize=8,
                   bbox_to_anchor=(1., -0.1))

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
        self.set_target()
        self.set_night_frames()
        self.set_target_list()
        self.set_plot()


if __name__ == "__main__":
    args = parse_args()
    Luke = Skywalker(args)
    Luke.main()
