#!/bin/python3

import os
import sys
import warnings
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body
from astropy.coordinates.errors import UnknownSiteException
from astropy.coordinates.name_resolve import NameResolveError
from astropy.time import Time
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import logging
from astroplan import Observer
from astroplan.plots import plot_sky
from timezonefinder import TimezoneFinder
import pytz

from .hover import TimeCursor, hover_is_possible
from .coords import parse_ra, parse_dec, parse_ra_column
from .plotdata import astro_night_mask


def parse_args(argv=None):
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
                        help="Right Ascension of the object, in degrees or \
                        hourangle, decimal or sexagesimal (e.g. '16:23:33.78' \
                        or '245.89'). The unit is auto-detected when \
                        possible; use --raunit to force it for an ambiguous \
                        value.")
    parser.add_argument("--dec", required=False, type=str,
                        help="Declination of the object in degrees. \
                        If using hexagesimal format, use --dec='DD:MM:SS'.")
    parser.add_argument("--raunit", type=str, default='auto',
                        help="Units of the RA parameter. Options: \
                        [auto, hour, deg]. 'auto' detects hourangle vs \
                        degrees per value and errors on a genuinely \
                        ambiguous one. Default is auto.")
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
    parser.add_argument("--no-hover", action="store_true",
                        help="Disable the interactive hover time cursor that \
                        shows the altitude and airmass of all targets at the \
                        time under the mouse pointer.")
    parser.add_argument("--savehtml", action="store_true",
                        help="Save an interactive HTML version of the figure, \
                        with hover enabled. Requires plotly. Default is False.")
    parser.add_argument("--htmlname", "-hn", type=str, required=False,
                        help="Name of the HTML file to save. \
                        Default is skywalker_<nightstarts>.html.")
    parser.add_argument("--htmljs", type=str, default="embed",
                        choices=["embed", "cdn", "directory"],
                        help="How to include plotly.js: embed (self-contained, \
                        works offline), cdn (small file, needs internet), or \
                        directory (shared plotly.min.js). Default is embed.")

    # web
    parser.add_argument("--web", action="store_true",
                        help="Serve an interactive web page with the plot \
                        and a target table below it, letting you include, \
                        remove or add targets live. Requires dash. Every \
                        browser tab connected shares the same target list. \
                        Default is False.")
    parser.add_argument("--web-host", type=str, default="127.0.0.1",
                        help="Interface to bind the web server to. \
                        Default is 127.0.0.1 (loopback only); a \
                        non-loopback host has no authentication in front \
                        of it, so prefer an SSH tunnel for remote access.")
    parser.add_argument("--web-port", type=int, default=8050,
                        help="TCP port for the web server. Default is 8050.")
    parser.add_argument("--web-open", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="Open the web UI in a new tab of the default \
                        browser once the server is listening, starting the \
                        browser if it is not running. On by default; use \
                        --no-web-open to disable.")
    parser.add_argument("--web-debug", action="store_true",
                        help="Enable the Dash/Werkzeug debugger. Refused \
                        unless --web-host is a loopback address, since the \
                        debugger allows remote code execution.")

    # logging
    parser.add_argument("--logfile", type=str, default='skywalker.log',
                        help="Log file to save the log messages.")
    parser.add_argument("--loglevel", type=str, default='INFO',
                        help="Log level. Options: [DEBUG, INFO, WARNING, ERROR, CRITICAL]. \
                        Default is INFO.")

    _argv = sys.argv[1:] if argv is None else argv
    if '-h' in _argv or '--help' in _argv:
        parser.print_help()
        sys.exit(0)
    try:
        args = parser.parse_args(argv)
    except SystemExit:
        raise
    except Exception:
        e = sys.exc_info()
        parser.error(f"Argument error: {e}. \
            If using hexagesimal DEC, use --dec='DD:MM:SS' to pass the argument.")
    if args.raunit not in ['auto', 'deg', 'hour']:
        parser.error(
            "Invalid value for --raunit. Options are: [auto, deg, hour].")
    if args.ra is not None:
        try:
            parse_ra(args.ra, raunit=args.raunit, context=" for --ra")
        except ValueError as e:
            parser.error(str(e))
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
        self.year_dates = None
        self.year_frame = None
        self.year_shape = None
        self.year_night_mask = None
        self.year_local_times = None
        self.year_of_frame = None

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
        self.make_hover = not args.no_hover
        self.savehtml = args.savehtml
        self.htmlname = args.htmlname
        self.htmljs = args.htmljs
        self.web = args.web
        self.web_host = args.web_host
        self.web_port = args.web_port
        self.web_open = args.web_open
        self.web_debug = args.web_debug
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.tracks = []
        self.time_cursor = None

        self.observer = None
        self.utcoffset = 0 * u.hour

        self.logfile = args.logfile
        self.loglevel = args.loglevel
        self.logger = logger(self.logfile, self.loglevel)

    def set_location(self, site=None, lat=None, lon=None, elev=None,
                     name=None):
        """Set the location of the observer based on the provided parameters.

        Parameters
        ----------
        site : str, optional
            Explicit site name from EarthLocation.get_site_names().
            Overridden by lat/lon when both of those are given. Falls
            back to self.site when None (the default).
        lat, lon : float, optional
            Explicit coordinates in degrees, given together. Checked
            first, so they take precedence over site when both are not
            None. Fall back to self.lat / self.lon when None (the
            default).
        elev : float, optional
            Elevation in meters, paired with lat/lon. Falls back to
            self.elev when None (the default).
        name : str, optional
            Display name for the site, used only when lat/lon are both
            given. Stripped of surrounding whitespace; None, or blank
            after stripping, falls back to the existing "{lat} {lon}"
            form. Ignored otherwise -- site, --sitefile and self.site
            already carry their own name.

        Any of these arguments bypasses --sitefile, which main()'s
        no-argument startup call still honours exactly as before this
        method gained parameters.
        """
        if lat is not None and lon is not None:
            _elev = elev if elev is not None else self.elev
            self.location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg,
                                          height=_elev * u.m)
            _name = (name or '').strip()
            self.sitename = _name if _name else f"{lat} {lon}"
        elif site is not None:
            try:
                self.location = EarthLocation.of_site(site)
            except UnknownSiteException:
                raise ValueError(
                    f"Site '{site}' not found in the database.")
            self.sitename = site
        elif self.sitefile is not None:
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
            except UnknownSiteException:
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
                    _time_str = f"{int(_time)}"
                    _time_str += f":{int((_time - int(_time)) * 60)}"
                    _time_str += ":00"
                    self.time = _time_str
            self.inithour = Time(self.nightstarts + "T" + self.time,
                                 format='isot').strftime('%H:%M:%S')
        tf = TimezoneFinder()
        timezone_str = tf.timezone_at(
            lng=self.location.lon.value, lat=self.location.lat.value)
        _timezone = pytz.timezone(timezone_str)
        self.utcoffset = (_timezone.utcoffset(Time(self.nightstarts + "T" +
                          self.inithour, format='isot').datetime).seconds / 3600 - 24) * u.hour
        self.obs_time = Time(self.nightstarts + "T" +
                             self.inithour, scale='utc', format='isot') - self.utcoffset

    def set_target(self, ra=None, dec=None):
        """Set the target object based on the provided parameters.
        If an object name is provided, it will try to resolve it to coordinates.
        If RA and DEC are provided, it will create a SkyCoord object.
        If a file is provided, it will load the coordinates from the file.

        Parameters:
        -----------
        ra : str or float, optional
            Right Ascension of the object, in degrees or hourangle, decimal
            or sexagesimal. The unit is auto-detected unless --raunit forces
            one; see coords.parse_ra().
        dec : str or float, optional
            Declination of the object in degrees.
        """
        if self.object:
            try:
                self.target = SkyCoord.from_name(self.object)
            except NameResolveError:
                raise ValueError(
                    f"Object '{self.object}' not found in the database.")
        elif (ra is not None) and (dec is not None):
            _ra_deg = parse_ra(ra, raunit=self.raunit, context=" for --ra")
            _dec_deg = parse_dec(dec)
            self.target = SkyCoord(ra=_ra_deg, dec=_dec_deg,
                                   unit=('deg', 'deg'))
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

    def set_year_frames(self, year=None, samples=145):
        """Build the shared time grid for the web UI's year-view panel.

        Builds one month-spaced night per calendar month of year -- the
        1st of January through the 1st of December -- each a 24 h window
        of samples epochs centred on that night's local midnight, all
        folded into a single vectorised AltAz frame. Requires
        set_location() and set_time() to have already run, since it
        relies on self.location, self.utcoffset and self.inithour.

        A call with year=None is a no-op once any grid has been built: it
        keeps whatever year is already on self.year_dates instead of
        forcing it back to the default, so year_max_altitudes()'s own
        no-argument call cannot silently undo a grid this method already
        built for a different year. Only an explicit year that differs
        from self.year_of_frame triggers a rebuild.

        Parameters:
        -----------
        year : int, optional
            Calendar year to build the grid for. None (the default) means
            "the year already built, if any, else the calendar year of
            self.nightstarts".
        samples : int, optional
            Number of epochs per night, spanning a 24 h window centred on
            local midnight. Default 145.
        """
        if self.year_frame is not None and (
                year is None or year == self.year_of_frame):
            return
        if year is not None:
            _year = year
        else:
            _year = pd.Timestamp(self.nightstarts).year

        # self.utcoffset was captured once for self.nightstarts (set_time(),
        # ~line 304) and is reused here for every month: nights half a year
        # away can end up centred up to an hour off of true local midnight
        # where DST applies. This is harmless because the window is 24 h
        # wide and we only take a maximum over it -- but this grid must
        # therefore never be repurposed for anything that displays clock
        # time to the user other than the coarse peak-time label computed in
        # year_max_altitudes().
        self.logger.debug(
            f"Building year frame grid: year={_year}, samples={samples}")

        _dates = [pd.Timestamp(year=_year, month=_m, day=1)
                 for _m in range(1, 13)]
        self.year_dates = [_d.date() for _d in _dates]
        self.year_of_frame = _year

        _midnight_list = []
        for _date in self.year_dates:
            _night_ends = (Time(f"{_date}T{self.inithour}", format='isot')
                          - self.utcoffset + .5 * u.day).strftime('%Y-%m-%d')
            _midnight = Time(f"{_night_ends}T00:00:00",
                             format='isot') - self.utcoffset
            _midnight_list.append(_midnight)
        _midnight_arr = Time(_midnight_list)

        _offsets = np.linspace(-12, 12, samples) * u.hour
        _flat_times = (_midnight_arr[:, None] + _offsets[None, :]).reshape(-1)

        self.year_frame = AltAz(obstime=_flat_times, location=self.location)
        _sun_alt_flat = get_body(
            'sun', _flat_times).transform_to(self.year_frame).alt.value
        _sun_alt_2d = _sun_alt_flat.reshape(12, samples)

        self.year_night_mask = astro_night_mask(_sun_alt_2d)
        self.year_local_times = (
            _flat_times + self.utcoffset).datetime.reshape(12, samples)
        self.year_shape = (12, samples)

    def year_max_altitudes(self, obj_coords):
        """Peak altitude, peak time and usable hours per month, per target.

        Altitude is maximised over astronomical night (Sun below -18 deg),
        which is tighter than the Sun-below-horizon (0 deg) cut used by
        compute_track() and by the web table's "Peak alt" column, so the two
        can legitimately disagree slightly for the same night -- that is not
        a bug.

        Parameters:
        -----------
        obj_coords : SkyCoord
            Coordinates of one target, or of N targets batched into a single
            SkyCoord by the caller (never looped in by this method).

        Returns:
        --------
        list of dict, one per target, in the same order as obj_coords:
            {'peak_alt': float ndarray (n_months,), NaN where not
                         observable,
             'peak_time': list[str] (n_months,), local 'HH:MM', '—' where
                          peak_alt is NaN,
             'hours_up': float ndarray (n_months,), hours above self.minalt
                         during astronomical night}
        """
        self.set_year_frames()

        _coords = obj_coords.reshape(1) if obj_coords.isscalar else obj_coords
        _n = _coords.size
        n_months, samples = self.year_shape

        _alt = _coords.reshape(-1, 1).transform_to(
            self.year_frame).alt.value.reshape(_n, n_months, samples)

        _masked = np.where(self.year_night_mask, _alt, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            _peak = np.nanmax(_masked, axis=2)
        _peak = np.where(_peak > 0., _peak, np.nan)

        _step_h = 24. / (samples - 1)
        _hours_up = (self.year_night_mask &
                    (_alt > self.minalt)).sum(axis=2) * _step_h

        _results = []
        for _t in range(_n):
            _peak_time = []
            for _m in range(n_months):
                if np.isnan(_peak[_t, _m]):
                    _peak_time.append('—')
                else:
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', RuntimeWarning)
                        _idx = np.nanargmax(_masked[_t, _m])
                    _peak_time.append(
                        self.year_local_times[_m, _idx].strftime('%H:%M'))
            _results.append({'peak_alt': _peak[_t],
                             'peak_time': _peak_time,
                             'hours_up': _hours_up[_t]})
        return _results

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
                _blinit_str = f"{int(_blinit)}"
                _blinit_str += f":{int((_blinit - int(_blinit)) * 60)}"
                _blinit_str += ":00"
                binit = _blinit_str
            _blockinit[i] = binit

        return _blockinit

    def parse_target_dataframe(self, df):
        """Validate and default-fill a target DataFrame, then parse it.

        Factored out of set_target_list()'s --file branch so the web UI's
        Load CSV button (webapp.TrackRegistry.replace_from_dataframe(),
        called from webapp._mutate()'s 'sw-csv-upload' branch) can apply
        exactly the same rules to an uploaded CSV without duplicating
        them: RA and DEC columns are required; NAME defaults to 'Obj1',
        'Obj2', ... ; BLINIT defaults to self.time and otherwise is
        normalised to 'HH:MM:SS' via check_blockinit_format() and
        anchored to self.nightstarts, exactly as an explicit BLINIT
        column always was; BLOCKTIME defaults to 0. RA is parsed with
        coords.parse_ra_column() (self.raunit) and DEC with
        coords.parse_dec() -- the same two calls set_target_list()
        always made.

        df is not mutated; a modified copy is returned. Raises
        ValueError, naming the missing column(s), when RA or DEC (or
        both) is absent. set_target_list() itself still checks this
        first, to keep its own file-naming message, so its call here
        never hits this branch; a caller with no file name of its own
        (the web upload path) relies on it directly.

        Returns a new DataFrame with NAME, RA (deg, float), DEC (deg,
        float), BLINIT ('HH:MM:SS' str) and BLOCKTIME (float, seconds)
        columns, row order preserved.
        """
        df = df.copy()
        _missing = [c for c in ('RA', 'DEC') if c not in df.columns]
        if _missing:
            raise ValueError(
                ' and '.join(_missing)
                + (' columns' if len(_missing) > 1 else ' column')
                + ' not found.')
        if 'NAME' not in df.columns:
            self.logger.warning(
                "NAME column not found. Using Obj incremented by 1 as "
                "names.")
            df['NAME'] = ["Obj" + str(i) for i in range(1, len(df) + 1)]
        if 'BLINIT' not in df.columns:
            self.logger.warning(
                "BLINIT column not found. Using time parameter as "
                "BLINIT.")
            df['BLINIT'] = [self.time] * len(df)
        else:
            df['BLINIT'] = self.check_blockinit_format(df['BLINIT'])
            df['BLINIT'] = df['BLINIT'].apply(
                lambda x: Time(self.nightstarts + "T" + x,
                              format='isot').strftime('%H:%M:%S'))
        if 'BLOCKTIME' not in df.columns:
            self.logger.warning(
                "BLOCKTIME column not found. Using 0 as BLOCKTIME.")
            df['BLOCKTIME'] = [0] * len(df)
        _ra_deg, _ra_warning = parse_ra_column(
            df['RA'], raunit=self.raunit, names=df['NAME'])
        if _ra_warning is not None:
            self.logger.warning(_ra_warning)
        df['RA'] = _ra_deg
        df['DEC'] = np.array([parse_dec(v) for v in df['DEC']])
        return df

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
            self.target_list = self.parse_target_dataframe(df)
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

    @staticmethod
    def set_empty_target_list():
        """An empty target list with the columns set_plot() expects."""
        return pd.DataFrame(columns=['NAME', 'RA', 'DEC', 'BLINIT',
                                     'BLOCKTIME'])

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

    def compute_track(self, name, ra, dec, blinit, blocktime=0., color=None):
        """Compute one target's track dict, without touching matplotlib.

        Parameters:
        -----------
        name : str
            Name of the target, as it should appear in the legend.
        ra, dec : float
            Coordinates in degrees, already resolved and unit-normalised.
        blinit : str
            Local start of the observing block, "HH:MM:SS".
        blocktime : float, optional
            Length of the observing block in seconds; 0 means no block.
        color : str, optional
            Colour to record for the renderers. None means "no colour of
            its own": set_plot() fills it in from the matplotlib artist and
            the HTML/hover renderers fall back to the shared palette.

        Returns the track dict, or None if the target is never above the
        horizon while the Sun is down -- the caller's cue to skip it.
        """
        obj_coords = SkyCoord(ra=ra, dec=dec, unit=('deg', 'deg'))

        block_starts = float(blinit.split(':')[0]) + \
            float(blinit.split(':')[1]) / 60. + \
            float(blinit.split(':')[2]) / 3600.
        if block_starts > 12.:
            block_starts -= 24

        myaltaz_overnight = obj_coords.transform_to(
            self.frame_time_overnight)

        mask = myaltaz_overnight.alt > 0 * u.deg
        mask &= self.sunaltaz_time_overnight.alt < 0 * u.deg
        if mask.sum() == 0:
            self.logger.warning(
                f"Object {name} is not observable at the given \
                    time and observatory.")
            return None

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

        if self.make_skychart:
            _chart_altaz = self.observer.altaz(observe_time, obj_coords)
            _chart_alt = _chart_altaz.alt.value
            _chart_az = _chart_altaz.az.value
        else:
            _chart_alt = _chart_az = None

        _block_ends = block_starts
        if blocktime > 0:
            _blocktime = float(blocktime) * u.s
            _blocktime = _blocktime.to(u.hour).value
            _block_ends = block_starts + _blocktime
        else:
            print("Block time is 0")

        moon_distance = self.moon.separation(obj_coords).value
        text_position = abs(self.delta_midnight.value - block_starts) == abs(
            self.delta_midnight.value - block_starts).min()
        if myaltaz_overnight.alt.value[text_position].size == 0:
            self.logger.warning(
                f"Could not find altitude for {name} at {block_starts}")
            altitude_position = 0.0
        elif myaltaz_overnight.alt.value[text_position].size > 1:
            self.logger.warning(
                f"Found more than one altitude for {name} at {block_starts}")
            altitude_position = myaltaz_overnight.alt.value[text_position].mean(
            )
        else:
            altitude_position = myaltaz_overnight.alt.value[text_position][0]

        moon_is_up = self.delta_midnight[self.moonaltaz_time_overnight.alt.value > 0].value
        if moon_is_up.size > 0 and (block_starts > moon_is_up.min()) and (block_starts < moon_is_up.max()):
            text_colour = 'magenta'
        else:
            text_colour = 'c'

        return {'name': str(name),
               'alt': myaltaz_overnight.alt.value,
               'az': myaltaz_overnight.az.value,
               'color': color,
               'is_moon': False,
               'has_block': blocktime > 0,
               'block_starts': block_starts,
               'block_ends': _block_ends,
               'moon_distance': moon_distance,
               'label_x': block_starts - 0.3,
               'label_y': altitude_position - 3,
               'label_color': text_colour,
               'chart_alt': _chart_alt,
               'chart_az': _chart_az,
               'chart_hours': hours_values if self.make_skychart else None,
               'coords': obj_coords,
               'observe_time': observe_time,
               'hours': hours_values}

    def set_plot(self):
        """Set the plot for the night observation."""
        self.tracks = []           # reset: set_plot may be called more than once
        if self.make_skychart:
            fig = self.fig = plt.figure(figsize=(16, 6))
            ax1 = self.ax1 = fig.add_subplot(121)
            ax3 = self.ax3 = fig.add_subplot(122, projection='polar')
        else:
            fig, ax1 = plt.subplots(figsize=(8, 6))
            self.fig, self.ax1 = fig, ax1
            ax3 = self.ax3 = None

        if self.target_list.index.size > 1:
            is_list = True
        else:
            is_list = False

        for index in self.target_list.index:
            myObjdf = self.target_list.loc[index]
            _track = self.compute_track(myObjdf['NAME'], myObjdf['RA'],
                                        myObjdf['DEC'], myObjdf['BLINIT'],
                                        blocktime=myObjdf['BLOCKTIME'])
            if _track is None:
                continue

            obj_coords = _track['coords']
            observe_time = _track['observe_time']
            hours_values = _track['hours']
            block_starts = _track['block_starts']
            block_ends = _track['block_ends']

            if is_list:
                p = ax1.plot(self.delta_midnight.value,
                             _track['alt'],
                             label=f"{myObjdf['NAME']}",
                             zorder=11)
                _mycolor = p[0].get_color()

                if myObjdf['BLOCKTIME'] > 0:
                    ax1.fill_between(self.delta_midnight.to('hr').value,
                                     np.zeros(
                        len(self.delta_midnight.value)),
                        _track['alt'],
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
                                 _track['alt'],
                                 c=_track['az'],
                                 label=myObjdf['NAME'],
                                 lw=0, s=8, cmap='viridis',
                                 zorder=11)
                plt.colorbar(sc, pad=0.1).set_label('Azimuth [deg]')
                _mycolor = None

                if myObjdf['BLOCKTIME'] > 0:
                    ax1.fill_between(self.delta_midnight.to('hr').value,
                                     np.zeros(
                        len(self.delta_midnight.value)),
                        _track['alt'],
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

            ax1.grid()

            ax1.text(_track['label_x'],
                     _track['label_y'],
                     "%i" % _track['moon_distance'],
                     fontsize=10, color=_track['label_color'], zorder=12)

            _track['color'] = _mycolor
            self.tracks.append(_track)

        ax1.plot(self.delta_midnight.to('hr').value,
                 self.moonaltaz_time_overnight.alt.value,
                 color='c', ls='--',
                 label='Moon: %i%%' % (
                     self.moon_brightness.value * 100),
                 zorder=10)
        ax1.fill_between(self.delta_midnight.to('hr').value, 0, 90,
                         (self.sunaltaz_time_overnight.alt < -0 *
                          u.deg) & (self.sunaltaz_time_overnight.alt > -6.3 * u.deg),
                         color='indigo', zorder=0, alpha=0.8)
        ax1.fill_between(self.delta_midnight.to('hr').value, 0, 90,
                         (self.sunaltaz_time_overnight.alt < -6 *
                          u.deg) & (self.sunaltaz_time_overnight.alt > -12.3 * u.deg),
                         color='indigo', zorder=1, alpha=0.9)
        ax1.fill_between(self.delta_midnight.to('hr').value, 0, 90,
                         (self.sunaltaz_time_overnight.alt < -12 *
                          u.deg) & (self.sunaltaz_time_overnight.alt > -18 * u.deg),
                         color='indigo', zorder=2, alpha=1)
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

        _moon_chart_alt = _moon_chart_az = _moon_chart_hours = None
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
                _moon_chart_altaz = self.observer.altaz(
                    moon_time, SkyCoord(ra=self.moon.ra, dec=self.moon.dec))
                _moon_chart_alt = _moon_chart_altaz.alt.value
                _moon_chart_az = _moon_chart_altaz.az.value
                _moon_chart_hours = moon_hours
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

        self.tracks.append({'name': 'Moon',
                            'alt': self.moonaltaz_time_overnight.alt.value,
                            'az': self.moonaltaz_time_overnight.az.value,
                            'color': 'c',
                            'is_moon': True,
                            'has_block': False,
                            'block_starts': 0.,
                            'block_ends': 0.,
                            'moon_distance': None,
                            'label_x': None,
                            'label_y': None,
                            'label_color': None,
                            'chart_alt': _moon_chart_alt,
                            'chart_az': _moon_chart_az,
                            'chart_hours': _moon_chart_hours})

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

        ax2 = self.ax2 = ax1.twinx()
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

        if self.savehtml:
            # Written before plt.show(), which blocks until the window closes.
            self.set_htmlplot()
        if self.make_hover:
            self.set_hover()
        plt.show()

    def set_hover(self):
        """Attach the interactive hover time cursor to the altitude plot."""
        if not self.tracks:
            self.logger.warning(
                "No observable targets to hover: cursor not enabled.")
            return
        if not hover_is_possible(self.fig):
            self.logger.info(
                f"Backend '{matplotlib.get_backend()}' is not interactive: \
                hover cursor not enabled.")
            return
        # Keep the reference: matplotlib holds callbacks weakly, so an
        # unreferenced cursor is garbage collected and hover silently stops.
        self.time_cursor = TimeCursor(self.fig, self.ax1, self.tracks,
                                      self.delta_midnight.value,
                                      ax2=self.ax2, ax3=self.ax3,
                                      utcoffset=self.utcoffset.value,
                                      minalt=self.minalt,
                                      logger=self.logger).connect()
        self.logger.info(
            f"Hover cursor enabled for {len(self.tracks)} track(s).")

    def _html_filename(self):
        """Name of the HTML file to save, mirroring the --figname convention."""
        if self.htmlname:
            name = self.htmlname
        else:
            name = f"skywalker_{self.nightstarts}.html"
        if not name.lower().endswith(('.html', '.htm')):
            name += '.html'
        return name

    def set_htmlplot(self):
        """Write an interactive HTML version of the figure, using plotly."""
        if not self.tracks:
            self.logger.warning(
                "No observable targets to plot: HTML figure not written.")
            return
        from . import htmlplot          # pure Python: no plotly needed to import it
        _local_times = (self.frame_time_overnight.obstime
                        + self.utcoffset).datetime
        try:
            _path = htmlplot.render(
                tracks=self.tracks,
                local_times=_local_times,
                delta_hours=self.delta_midnight.value,
                sun_alt=self.sunaltaz_time_overnight.alt.value,
                moon_alt=self.moonaltaz_time_overnight.alt.value,
                moon_brightness=self.moon_brightness.value,
                minalt=self.minalt,
                utcoffset_h=int(self.utcoffset.value),
                sitename=self.sitename,
                nightstarts=self.nightstarts,
                make_skychart=self.make_skychart,
                filename=self._html_filename(),
                include_js=self.htmljs,
                logger=self.logger)
        except ImportError as e:
            self.logger.error(str(e))
            return
        self.logger.info(f"Interactive figure saved as {_path}")
        print(f"Interactive figure saved as {_path}")

    def set_webapp(self):
        """Serve the interactive plot and target table over HTTP, using dash."""
        matplotlib.use('Agg')      # web mode never opens a GUI figure
        from . import webapp
        webapp.run_webapp(self, host=self.web_host, port=self.web_port,
                          debug=self.web_debug, open_browser=self.web_open)

    def main(self):
        if self.web:
            self.set_location()
            self.set_observer()
            self.set_time()
            self.set_night_frames()
            if (self.object or self.file
                    or (self.ra is not None and self.dec is not None)):
                self.set_target(ra=self.ra, dec=self.dec)
                self.set_target_list()
            else:
                self.target_list = Skywalker.set_empty_target_list()
            self.set_webapp()
            return

        self.set_location()
        self.set_observer()
        self.set_time()
        self.set_target(ra=self.ra, dec=self.dec)
        self.set_night_frames()
        self.set_target_list()
        self.set_plot()


def main():
    args = parse_args()
    if args.sites:
        print("Available sites:")
        for site in EarthLocation.get_site_names():
            print(site)
    else:
        Luke = Skywalker(args)
        Luke.main()


if __name__ == "__main__":
    main()
