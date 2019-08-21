#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon, get_sun
import sys, os
import datetime
from astropy.io import ascii
import pytz
from tzwhere import tzwhere
import warnings
from myastroplan.plots import plot_sky
from myastroplan import Observer

lastversion = '0.4'
moddate = '2019-08-16'

initext = """
    ===================================================================
                    SkyWalker - v%s - %s
             This version is not yet completely debugged
             In case of crashes, please send the log to:
                Herpich F. R. fabiorafaelh@gmail.com
    ===================================================================
    """ % (lastversion, moddate)
print(initext)

def print_help():

    text_help = """

    SkyWalker - the skycalculator fodastico!

    skywalker.py - v%s - %s
    this version was recently modified and not yet completely debugged.
    In case of crashes, please send the log to:
    Herpich F. R. fabiorafaelh@gmail.com

    Usage: skywalker.py [--option] [entry]

    Calculates the target(s) track(s) during the night.

    --help, -h:             print this help

    --nightstarts, -ns      [YY-MM-DD] date of the starting night. Default is today
                            Other options to replace with:
    --today                 shows track for today
    --tomorrow              shows track for tomorrow
    --yesterday             shows track from yesterday

    --at                    [hh:mm:ss] shows moon brigthness and angular distance
                            to object at the specified hour. This value is also used
                            to start the observation block for an especified object
                            when defined. Must be Local Time
    --astrotime             show the total number of astronomical hours of the night

    --object                [string] to entry a known object by name. It will look at the
                            internal database for it
    --name (optional)       [string] a user defined name for the object
    --names                 [string1,string2,...] if you want to enter a list of
                            objects through command line, this option is mandatory.
                            The names should be separated by comma
    --ra                    [value|value1,value2,...] if list mode --names option
                            is selected, the entry for --ra should be a sequence
                            of values separated by comma with the same quantity
                            as for --names. Otherwise, a single value should be given
    --dec                   [value|value1,value2,...] same as --ra
    --file, -f              [string] should be a list containing the names of the
                            targets (optional), its RA (mandatory) and DEC (mandatory)
                            in this order. The file should be in ascii mode with
                            the entrys separated by space or comma.
                            The file can contain as many columns the user wants, but
                            in such a case, the existence of the columns named
                            NAME, RA and DEC (upercase) is mandatory.
                            The file containing the list of fields can contain two
                            other columns (they are not mandatory). The columns should
                            be named BLINIT (initial time in hours) and the other
                            BLOCKTIME (in seconds). The first will set the initial
                            time to start tracking the given object. Is this value
                            that will be used to calculate the distance to the moon.
                            The second sets the length of the observation block. If
                            BLOCKTIME is not given, it will set the same length for all
                            targets in the list using the --blocktime entry. If neither
                            --blocktime and BLOCKTIME are given, the initial time
                            If both --blocktime and BLOCKTIME are given, --blocktime
                            paramter is ignored and only BLOCKTIME is used.
                            BLINIT will be used only for moon distance calculation.
    --pid                   [string] if file is provided, --pid can be used to select
                            only part of the objects in it. It is mandatory the file
                            to have the column PID. If it hasn't, the argument --pid
                            is ignored. When provided, only the file that match the
                            argument given by --pid will be shown.

    --observatory           [string] the name of a known observatory. It will look
                            for it at the database (an internet connection is required).
                            To see the list of observatory, just enter this option
                            with no other argument. Default is Cerro Tololo.
    --lat                   [value] latitude for a user defined observatory
    --lon                   [value] longitude for a user defined observatory
    --height                [value] height for a user defined observatory in
                            meters. Default is 1000m
    --sitefile              [string] a csv file containing all the information
                            about the site. Can be used instead of using --observatory,
                            --lon and height. The columns must be named as NAME,LON,LAT,HEIGHT
    --minalt                [float] altitude limit of the telescope intended to be used
                            for the observations

    --skychart              Define the track of the target in the sky starting at
                            the time defined by --at

    --savechart             [string] if the user want to save the chart plan. The name
                            here will be used to save the plan into the direcotyr figs/

    --blocktime             [float] time that an entire observation block should
                            take (in seconds). Only shown for single objects.
    --noshow                do not show plt.figure
    """ % (lastversion, moddate)
    print(text_help)
    sys.exit(0)


def skychart(ax, mysite, mycoords, time, observe_time, sunaltazs_time_overnight,
             obj_style={'color': 'b'}, hours_value=False):
    """
    draws the map of the sky for the night
    """

    observer = Observer(longitude=mysite.lon, latitude=mysite.lat, elevation=mysite.height)

    plot_sky(mycoords, observer, observe_time, ax=ax, style_kwargs=obj_style, hours_value=hours_value)

    return ax

def standard_hour(inithour, night_starts, night_ends):
    """
    define the proper hour for the night starting
    """

    if len(inithour.split(':')) == 1:
        inithour += ':00:00'
        print('INFO - estimating distance to the moon at %s' % inithour)
    elif len(inithour.split(':')) == 2:
        inithour += ':00'
        print('INFO - estimating distance to the moon at %s' % inithour)
    elif (len(inithour.split(':')) < 1) | (len(inithour.split(':')) > 3):
        raise SyntaxError('wrong hour format. Use hh:mm:ss')
    if (np.float(inithour.split(':')[0]) >= 0.) & (np.float(inithour.split(':')[0]) < 12):
        night_starts = night_ends
    else:
        night_starts = night_starts
    test_value = 0
    for hval in inithour.split(':'):
        test_value += np.int(hval)
    if test_value == 0:
        inithour = '00:00:00.1'

    return inithour, night_starts

####################################

if len(sys.argv) < 2:
    print_help()
elif ('--help' in sys.argv) | ('-h' in sys.argv):
    print_help()
elif (len(sys.argv) == 2) & ('--observatory' in sys.argv):
    for site in EarthLocation.get_site_names(): print(site)
    sys.exit(0)

# define the starting night
if ('--nightstarts' in sys.argv) | ('-ns' in sys.argv):
    if ('--nightstarts' in sys.argv):
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--nightstarts']
    else:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '-ns']
    night_starts = np.array(sys.argv)[index + 1].item()
    night_beg = night_starts
    if len(night_starts.split('-')) < 3:
        raise IOError('wrong format for date')
    else:
        startsin = Time('%sT23:59:59' % night_starts)
        ns = startsin.datetime
        night_starts = ns.strftime("%Y-%m-%d")
elif '--tomorrow' in sys.argv:
    ns = datetime.datetime.now() + datetime.timedelta(days=1)
    night_starts = ns.strftime("%Y-%m-%d")
elif '--yesterday' in sys.argv:
    ns = datetime.datetime.now() - datetime.timedelta(days=1)
    night_starts = ns.strftime("%Y-%m-%d")
elif '--today' in sys.argv:
    ns = datetime.datetime.now()
    night_starts = ns.strftime("%Y-%m-%d")
else:
    ns = datetime.datetime.now()
    night_starts = ns.strftime("%Y-%m-%d")
ne = ns + datetime.timedelta(days=1)
night_ends = ne.strftime("%Y-%m-%d")

# build the necessary requisites to draw the observation block
if '--blocktime' in sys.argv:
    if '--at' in sys.argv:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--blocktime']
        _blocktime = np.array(sys.argv)[index + 1].item()
        try:
            blocktime = np.float(_blocktime)
            bl_time = True
        except:
            warnings.warn('WARNING - wrong time format for blocktime (should be seconds)')
            bl_time = False
    else:
        warnings.warn('WARNING - initial time for block not defined. Skipping')
        bl_time = False
else:
    bl_time = False

# deal with the input list
if ('-f' in sys.argv) | ('--file' in sys.argv):
    isList = True
    isKnown=False
    # check if it is all ok with the list
    if ('-f' in sys.argv) & ('--file' in sys.argv):
        print_help()
        raise IOError('you should chose either -f or --file option')
    elif ('-f' in sys.argv) & ('--file' not in sys.argv):
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '-f']
    else:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--file']
    myListObj = []
    mytablename = np.array(sys.argv)[index + 1].item()
    f = ascii.read(mytablename)
    if len(f.keys()) < 2:
        raise IOError('wrong format for list. Needs at least two columns comma or space separated')
    elif len(f.keys()) == 2:
        noname = True
        racol = f[f.keys()[0]]
        deccol = f[f.keys()[1]]
    elif len(f.keys()) == 3:
        noname = False
        namecol = f[f.keys()[0]]
        racol = f[f.keys()[1]]
        deccol = f[f.keys()[2]]
    elif len(f.keys()) > 3:
        warn_text = 'for file with more than 3 columns is mandatory to have the'
        warn_text += ' columns NAME, RA and DEC'
        warnings.warn(warn_text)
        if 'NAME' not in f.keys():
            raise IOError('column NAME not found')
        elif 'RA' not in f.keys():
            raise IOError('column RA not found')
        elif 'DEC' not in f.keys():
            raise IOError('column DEC not found')
        noname = False
        namecol = f['NAME']
        racol = f['RA']
        deccol = f['DEC']
    # look for some optional cols
    if 'BLINIT' in f.keys():
        blinit = True
    else:
        blinit = False
    if 'BLOCKTIME' in f.keys():
        bllenght = True
    else:
        bllenght = False
    # prepare to plot a list of specific objects (PIDs) inside a bigger list
    if '--pid' in sys.argv:
        if 'PID' not in f.keys():
            warnings.warn('no PID column found. Proceeding for all objects within the file')
            maskpid = np.full(racol.size, True)
        else:
            index = np.arange(len(sys.argv))[np.array(sys.argv) == '--pid']
            pid = np.array(sys.argv)[index + 1].item()
            if '--todo' in sys.argv:
                maskpid = (f['PID'] == pid) & ((f['STATUS'] == 0) | (f['STATUS'] == 3))
            else:
                maskpid = f['PID'] == pid
            print('INFO - calculating for', pid)
    elif '--name' in sys.argv:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--name']
        pid = np.array(sys.argv)[index + 1].item()
        maskpid = f['NAME'] == pid
        # draw the observation block for a list even with no BLOCKTIME column
        # present in list
        if '--blocktime' in sys.argv:
            if '--at' in sys.argv:
                index = np.arange(len(sys.argv))[np.array(sys.argv) == '--at']
                inithour = np.array(sys.argv)[index + 1].item()
                index = np.arange(len(sys.argv))[np.array(sys.argv) == '--blocktime']
                _blocktime = np.array(sys.argv)[index + 1].item()
                try:
                    blocktime = np.float(_blocktime)
                    bl_time = True
                except:
                    warnings.warn('WARNING - wrong time format for blocktime (should be seconds)')
                    bl_time = False
            else:
                warnings.warn('WARNING - initial time for block not defined. Skipping')
                bl_time = False
    else:
        maskpid = np.full(racol.size, True)

    # build the list of objects to be ploted
    for j in range(racol[maskpid].size):
        myListObj.append([])
        if noname:
            myListObj[j].append('obj' + str(j + 1))
        else:
            myListObj[j].append(namecol[maskpid][j])
        if (len(str(racol[maskpid][j]).split(':')) > 1) | (len(str(racol[maskpid][j]).split('h')) > 1):
            raunit = u.hourangle
        else:
            raunit = u.deg
        objcoords = SkyCoord(racol[maskpid][j], deccol[maskpid][j], unit=(raunit, u.deg))
        myListObj[j].append(objcoords.ra.value)
        myListObj[j].append(objcoords.dec.value)
        if blinit:
            myListObj[j].append(f['BLINIT'][maskpid][j])
            make_blinit = True
        elif bl_time:
            myListObj[j].append(inithour)
            make_blinit = True
            blinit = True
        else:
            make_blinit = False
        if bllenght:
            myListObj[j].append(f['BLOCKTIME'][maskpid][j])
        elif (not bllenght) & bl_time:
            myListObj[j].append(blocktime)
        else:
            pass
elif '--names' in sys.argv:# if list is not the given, but using sequence of names
    warnings.warn('using names, you should provide the names and positions separated by coma')
    isList = True
    isKnown=False
    index = np.arange(len(sys.argv))[np.array(sys.argv) == '--names']
    myobjects = np.array(sys.argv)[index + 1].item().split(',')
    index = np.arange(len(sys.argv))[np.array(sys.argv) == '--ra']
    myras = np.array(sys.argv)[index + 1].item().split(',')
    index = np.arange(len(sys.argv))[np.array(sys.argv) == '--dec']
    mydecs = np.array(sys.argv)[index + 1].item().split(',')
    myListObj = []
    for j in range(len(myobjects)):
        myListObj.append([])
        myListObj[j].append(myobjects[j])
        if len(myras[j].split(':')) > 1:
            raunit = u.hourangle
        else:
            raunit = u.deg
        objcoords = SkyCoord(myras[j], mydecs[j], unit=(raunit, u.deg))
        myListObj[j].append(objcoords.ra.value)
        myListObj[j].append(objcoords.dec.value)
    blinit = False
elif '--object' in sys.argv:# if single object is provided
    isKnown = True
    isList = False
    index = np.arange(len(sys.argv))[np.array(sys.argv) == '--object']
    objname = np.array(sys.argv)[index + 1].item()
    myListObj = [[objname]]
    blinit = False
else: # otherwise it is mandatory entering ra and dec
    isKnown=False
    isList = False
    if (('--ra' not in sys.argv) | ('--dec' not in sys.argv)) & ('--coords' not in sys.argv):
        print_help()
        raise SyntaxError('check your options')
    elif ('--name' not in sys.argv):
        indexes = np.arange(len(sys.argv))
        objname = 'obj01'
        if '--coords' in sys.argv:
            if len(np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 1].item().split(',')) > 1:
                objra = np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 1].item().split(',')[0]
                objdec = np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 1].item().split(',')[1]
            elif len(np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 1].item().split(' ')) > 1:
                objra = np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 1].item().split(' ')[0]
                objdec = np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 1].item().split(' ')[1]
            else:
                objra = np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 1].item()
                objdec = np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 2].item()
        else:
            objra = np.array(sys.argv)[indexes[np.array(sys.argv) == '--ra'] + 1].item()
            objdec = np.array(sys.argv)[indexes[np.array(sys.argv) == '--dec'] + 1].item()
        if len(objra.split(':')) > 1:
            raunit = u.hourangle
        else:
            raunit = u.deg
        objcoords = SkyCoord(objra, objdec, unit=(raunit, u.deg))
    else:
        indexes = np.arange(len(sys.argv))
        objname = np.array(sys.argv)[indexes[np.array(sys.argv) == '--name'] + 1].item()
        if '--coords' in sys.argv:
            objra = np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 1].item()
            objdec = np.array(sys.argv)[indexes[np.array(sys.argv) == '--coords'] + 2].item()
        else:
            objra = np.array(sys.argv)[indexes[np.array(sys.argv) == '--ra'] + 1].item()
            objdec = np.array(sys.argv)[indexes[np.array(sys.argv) == '--dec'] + 1].item()
        if len(objra.split(':')) > 1:
            raunit = u.hourangle
        else:
            raunit = u.deg
        objcoords = SkyCoord(objra, objdec, unit=(raunit, u.deg))
    myListObj = [[objname, objcoords.ra.value, objcoords.dec.value]]
    blinit = False

# check for site
if '--observatory' in sys.argv:# if site is provided
    index = np.arange(len(sys.argv))[np.array(sys.argv) == '--observatory']
    try:
        sitename = np.array(sys.argv)[index + 1].item()
    except:
        print('please choose an observatory among the available at the')
        print('following list. Otherwise, provide positional data for location')
        for site in EarthLocation.get_site_names(): print(site)
        sys.exit(0)
    if (sitename not in EarthLocation.get_site_names()) & (('--lat' not in sys.argv) & ('--lon' not in sys.argv)):
        print('please choose an observatory among the available at the')
        print('following list. Otherwise, provide positional data for location')
        for site in EarthLocation.get_site_names(): print(site)
        sys.exit(0)
    elif (sitename not in EarthLocation.get_site_names()) & ('--lat' in sys.argv) & ('--lon' in sys.argv):
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--lat']
        lat = np.array(sys.argv)[index + 1].item()
        if lat.split(':') > 1:
            lapos = lat.split(':')
            if lapos[0].split('-') > 1:
                multiplier = -1.
            else:
                multiplier = 1.
            _mylat = 0.
            divisor = 1.
            for val in lapos:
                _mylat += abs(np.float(val)) / divisor
                divisor *= 60.
            mylat = multiplier * _mylat
        else:
            mylat = np.float(lat)
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--lon']
        lon = np.array(sys.argv)[index + 1].item()
        if lon.split(':') > 1:
            lopos = lon.split(':')
            if lopos[0].split('-') > 1:
                multiplier = -1.
            else:
                multiplier = 1.
            _mylon = 0.
            divisor = 1.
            for val in lopos:
                _mylon += abs(np.float(val)) / divisor
                divisor *= 60.
            mylon = multiplier * _mylon
        if '--height' in sys.argv:
            index = np.arange(len(sys.argv))[np.array(sys.argv) == '--height']
            height = np.float(np.array(sys.argv)[index + 1].item())
        else:
            warnings.warn('height not provided. Using default 1000 m')
            height = 1000.
        mysite = EarthLocation(lat=mylat*u.deg, lon=mylon*u.deg, height=height*u.m)
        if sitename.split('-')[0] == '':
            sitename = 'LAT%.1f_LON%.1f' % (mylat, mylon)
    elif sitename in EarthLocation.get_site_names():
        mysite = EarthLocation.of_site(sitename)
    tzwhere = tzwhere.tzwhere()
    timezone_str = tzwhere.tzNameAt(mysite.lat.value, mysite.lon.value)
    timezone = pytz.timezone(timezone_str)
    dt = datetime.datetime(int(night_starts.split('-')[0]),
                           int(night_starts.split('-')[1]),
                           int(night_starts.split('-')[2]),
                           0, 0, 0, 0)
    # only UTC is implemented
    utcoffset = (timezone.utcoffset(dt).seconds/3600 - 24)*u.hour
    if (utcoffset.value < 0) & (utcoffset.value >= -12):
        xlabel = 'Local Time (UTC %i)' % utcoffset.value
    else:
        xlabel = 'Local Time (UTC + %i)' % (utcoffset.value + 24)
    tlabel = 'LT'
    #utcoffset = 0*u.hour
    #xlabel = 'UTC'
    #tlabel = 'UTC'
    if '--minalt' in sys.argv:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--minalt']
        minalt = np.float(np.array(sys.argv)[index + 1].item())
    else:
        minalt = 0
elif '--sitefile' in sys.argv:# if file containing site info is given
    index = np.arange(len(sys.argv))[np.array(sys.argv) == '--sitefile']
    sitefilename = np.array(sys.argv)[index + 1].item()
    sitefile = ascii.read(sitefilename)
    sitename = sitefile['NAME'][0]
    mysite = EarthLocation(lat=sitefile['LAT'][0], lon=sitefile['LON'][0], height=sitefile['HEIGHT'][0])
    #utcoffset = 0*u.hour
    #xlabel = 'UTC'
    #tlabel = 'UTC'
    if '--minalt' in sys.argv:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--minalt']
        minalt = np.float(np.array(sys.argv)[index + 1].item())
    else:
        minalt = 0
    tzwhere = tzwhere.tzwhere()
    timezone_str = tzwhere.tzNameAt(mysite.lat.value, mysite.lon.value)
    timezone = pytz.timezone(timezone_str)
    dt = datetime.datetime(int(night_starts.split('-')[0]),
                           int(night_starts.split('-')[1]),
                           int(night_starts.split('-')[2]),
                           0, 0, 0, 0)
    # only Local Time is implemented
    utcoffset = (timezone.utcoffset(dt).seconds/3600 - 24)*u.hour
    if (utcoffset.value < 0) & (utcoffset.value >= -12):
        xlabel = 'Local Time (UTC %i)' % utcoffset.value
    else:
        xlabel = 'Local Time (UTC + %i)' % (utcoffset.value + 24)
    tlabel = 'LT'
else:# the default is the Cerro Tololo Observatory, with some specs used by the S-PLUS-T80S team
    warnings.warn('Observatory not identified. Using default T80S at Cerro Tololo')
    sitename = 'T80-South'
    mysite = EarthLocation(lat=-30.2*u.deg, lon=-70.8*u.deg, height=2200*u.m)
    #utcoffset = 0*u.hour
    #xlabel = 'UTC'
    #tlabel = 'UTC'
    minalt = 30
    tzwhere = tzwhere.tzwhere()
    timezone_str = tzwhere.tzNameAt(mysite.lat.value, mysite.lon.value)
    timezone = pytz.timezone(timezone_str)
    dt = datetime.datetime(int(night_starts.split('-')[0]),
                           int(night_starts.split('-')[1]),
                           int(night_starts.split('-')[2]),
                           0, 0, 0, 0)
    # only Local Time is implemented
    utcoffset = (timezone.utcoffset(dt).seconds/3600 - 24)*u.hour
    if (utcoffset.value < 0) & (utcoffset.value >= -12):
        xlabel = 'Local Time (UTC %i)' % utcoffset.value
    else:
        xlabel = 'Local Time (UTC + %i)' % (utcoffset.value + 24)
    tlabel = 'LT'

night_for_chart = night_starts
# define the time to calculate the distance to the moon and to start the
# observing block
if '--at' in sys.argv:
    index = np.arange(len(sys.argv))[np.array(sys.argv) == '--at']
    inithour = np.array(sys.argv)[index + 1].item()
    #if inithour >= '24':
    #    raise IOError('wrong initial hour %s' % inithour)
    inithour, night_starts = standard_hour(inithour, night_starts, night_ends)
else:
    inithour = '23:59:59'

time = Time('%s %s' % (night_starts, inithour)) - utcoffset
titlenight = 'Conditions at ' + (time + utcoffset).value[:-4] + ' LT for ' + sitename
midnight = Time('%s 00:00:00' % night_ends) - utcoffset
delta_midnight = np.linspace(-12, 12, 500)*u.hour
frame_tonight = AltAz(obstime=midnight+delta_midnight,
                          location=mysite)
dec_inhour = np.float(inithour.split(':')[0])
dec_inhour += np.float(inithour.split(':')[1]) / 60.
dec_inhour += np.float(inithour.split(':')[2]) / 3600.
if dec_inhour > 12.:
    dec_inhour = dec_inhour - 24.

# coords frames
times_time_overnight = midnight + delta_midnight
frame_time_overnight = AltAz(obstime=times_time_overnight, location=mysite)
sunaltazs_time_overnight = get_sun(times_time_overnight).transform_to(frame_time_overnight)

moon_time_overnight = get_moon(times_time_overnight)
moonaltazs_time_overnight = moon_time_overnight.transform_to(frame_time_overnight)

# moon phase
sun_pos = get_sun(time+utcoffset)
moon_pos = get_moon(time+utcoffset)
elongation = sun_pos.separation(moon_pos)
moon_phase = np.arctan2(sun_pos.distance*np.sin(elongation),
                        moon_pos.distance - sun_pos.distance*np.cos(elongation))
moon_brightness = (1. + np.cos(moon_phase))/2.0
# transparency to the background of the plots
if moon_brightness.value < .1:
    _moon_brightness = .1
else:
    _moon_brightness = moon_brightness.value
mb = int(moon_brightness*100)

# start the figure
if '--skychart' in sys.argv:
    fig = plt.figure(figsize=(16,6))
    ax1 = fig.add_subplot(1,2,1)
    ax3 = fig.add_subplot(1,2,2, projection='polar')
else:
    fig, ax1 = plt.subplots()

# start building the figure
for myObj in myListObj:
    # define the object coords, time, alt, az
    if isKnown:
        mycoords = SkyCoord.from_name(objname)
    else:
        mycoords = SkyCoord(ra=myObj[1]*u.deg, dec=myObj[2]*u.deg)

    if blinit:
        objhour, night_starts = standard_hour(myObj[3], night_starts, night_ends)
    else:
        objhour = inithour
    objtime = Time('%s %s' % (night_starts, objhour)) - utcoffset
    obj_dec_inhour =  np.float(objhour.split(':')[0])
    obj_dec_inhour += np.float(objhour.split(':')[1]) / 60.
    obj_dec_inhour += np.float(objhour.split(':')[2]) / 3600.
    if obj_dec_inhour > 12.:
        obj_dec_inhour = obj_dec_inhour - 24.

    inivalue = obj_dec_inhour
    endvalue = delta_midnight[sunaltazs_time_overnight.alt < -18*u.deg].max().value
    #observe_time = Time('%s 23:59:59.9' % night_for_chart) + np.arange(inivalue, endvalue, 1)*u.hour
    observe_time = time + np.arange(inivalue, endvalue, 1)*u.hour

    myaltaz = mycoords.transform_to(AltAz(obstime=objtime, location=mysite))

    myaltazs_tonight = mycoords.transform_to(frame_tonight)

    myaltazs_time_overnight = mycoords.transform_to(frame_time_overnight)

    if isList:# if list
        p = ax1.plot(delta_midnight, myaltazs_time_overnight.alt, label=myObj[0])
        titlenight = 'Conditions for night starting at ' + (time + utcoffset).value[:10] + ' LT for ' + sitename
        if '--skychart' in sys.argv:
            skychart(ax3, mysite, mycoords, time, observe_time,
                     sunaltazs_time_overnight, obj_style={'color': p[0].get_color(),
                                                          'marker': '*',
                                                          'label': myObj[0]}, hours_value=np.arange(inivalue, endvalue, 1))
        if (bl_time or bllenght) & blinit:
            try:
                obj_bltime = np.float(myObj[4]) / 3600.
                ax1.fill_between(delta_midnight,
                                 np.zeros(len(delta_midnight)),
                                 myaltazs_time_overnight.alt,
                                 (delta_midnight.value >= obj_dec_inhour) & (delta_midnight.value <= (obj_dec_inhour + obj_bltime)),
                                 color=p[0].get_color())
            except:
                warnings.warn('wrong time format for blocktime (should be seconds)')
    else:# if single object
        sc = ax1.scatter(delta_midnight, myaltazs_time_overnight.alt,
                    c=myaltazs_time_overnight.az, label=myObj[0], lw=0, s=8,
                    cmap='viridis')
        if bl_time:
            ax1.fill_between(delta_midnight,
                             np.zeros(len(delta_midnight)),
                             myaltazs_time_overnight.alt,
                             (delta_midnight.value >= dec_inhour) & (delta_midnight.value <= (dec_inhour + blocktime / 3600.)),
                             color='y')
        if '--skychart' in sys.argv:
            skychart(ax3, mysite, mycoords, time+utcoffset, observe_time,
                     sunaltazs_time_overnight, obj_style={'cmap': 'viridis',
                                                          'marker': '*',
                                                          'c': np.arange(inivalue, endvalue, 1),
                                                          'label': myObj[0]}, hours_value=np.arange(inivalue, endvalue, 1))

    text_pos = abs(delta_midnight.value - obj_dec_inhour) == abs(delta_midnight.value - obj_dec_inhour).min()
    ax1.text(obj_dec_inhour - 0.3,
             myaltazs_time_overnight.alt.value[text_pos] - 3.0,
             '%i' % moon_pos.separation(mycoords).value, fontsize=10, color='c')

# plot moon and the rest of the paramentrs
ax1.plot(delta_midnight, moonaltazs_time_overnight.alt, color='violet', ls='--',
         label=r'Moon: %i%%' % mb)
ax1.fill_between(delta_midnight.to('hr').value, 0, 90,
                 (sunaltazs_time_overnight.alt < -0*u.deg) & (sunaltazs_time_overnight.alt > -18*u.deg),
                 color='indigo', zorder=0)
ax1.fill_between(delta_midnight.to('hr').value, 0, 90,
                 (moonaltazs_time_overnight.alt <= 0*u.deg) & (sunaltazs_time_overnight.alt < -18*u.deg), color='k', zorder=0)
ax1.fill_between(delta_midnight.to('hr').value, 0, 90,
                 (moonaltazs_time_overnight.alt > 0*u.deg) & (sunaltazs_time_overnight.alt < -18*u.deg),
                 color='midnightblue', alpha=1.1 - _moon_brightness, zorder=0)

if '--skychart' in sys.argv:
    skychart(ax3, mysite, SkyCoord(ra=moon_pos.ra, dec=moon_pos.dec), time, observe_time,
             sunaltazs_time_overnight, obj_style={'color': 'cyan',
                                                  'marker': 'o',
                                                  'label': 'moon: %i%%' % (moon_brightness*100)},
             hours_value=np.arange(inivalue, endvalue, 1))
    circle = plt.Circle((0., 0.), 90, transform=ax3.transData._b,
                       color="red", alpha=0.7, zorder=0)
    ax3.add_artist(circle)
    circle = plt.Circle((0., 0.), 90-minalt, transform=ax3.transData._b,
                       color="white", alpha=1., zorder=0)
    ax3.add_artist(circle)
    circle = plt.Circle((0., 0.), 90-minalt, transform=ax3.transData._b,
                       color="black", alpha=1.1 - _moon_brightness, zorder=0)
    ax3.add_artist(circle)
if isList:
    if '--astrotime' in sys.argv:
        unt = delta_midnight.value[sunaltazs_time_overnight.alt <= -18*u.deg].max()
        unt -= delta_midnight.value[sunaltazs_time_overnight.alt <= -18*u.deg].min()
        if len(myListObj[0]) <= 20:
            if '--skychart' in sys.argv:
                leg = plt.legend(loc='lower left', title='Astronomical time: %.2fh' % unt,
                ncol=2, fancybox=True, fontsize=11, bbox_to_anchor=(.8, 0.), scatterpoints=1)
            else:
                leg = plt.legend(loc='lower left', title='Astronomical time: %.2fh' % unt,
                ncol=4, fancybox=True, fontsize=11)
            leg.get_frame().set_alpha(0.9)
            leg.get_frame().set_edgecolor('white')
    else:
        if len(myListObj) <= 20:
            if '--skychart' in sys.argv:
                leg = plt.legend(loc='lower left', ncol=2, fancybox=True,
                                 fontsize=11, bbox_to_anchor=(.8, -0.2), scatterpoints=1)
            else:
                leg = plt.legend(loc='lower left', ncol=4, fancybox=True, fontsize=11)
            leg.get_frame().set_alpha(0.9)
            leg.get_frame().set_edgecolor('white')
else:
    if '--astrotime' in sys.argv:
        unt = delta_midnight.value[sunaltazs_time_overnight.alt <= -18*u.deg].max()
        unt -= delta_midnight.value[sunaltazs_time_overnight.alt <= -18*u.deg].min()
        if '--skychart' in sys.argv:
            leg = plt.legend(loc='center left', title='Astronomical time: %.2fh' % unt,
            fancybox=True, scatterpoints=1, bbox_to_anchor=(.8, 0.))
        else:
            leg = plt.legend(loc='lower left', title='Astronomical time: %.2fh' % unt,
            fancybox=True, scatterpoints=1)
    else:
        if '--skychart' in sys.argv:
            leg = plt.legend(loc='center left', scatterpoints=1, bbox_to_anchor=(1., 0.), fancybox=True)
        else:
            leg = plt.legend(loc='lower left', fancybox=True, scatterpoints=1)
    leg.get_frame().set_alpha(0.9)
    leg.get_frame().set_edgecolor('white')
minx = delta_midnight.value[sunaltazs_time_overnight.alt < -0*u.deg].min() - 1
maxx = delta_midnight.value[sunaltazs_time_overnight.alt < -0*u.deg].max() + 1
if minalt > 1:
    ax1.plot([minx, maxx], [minalt, minalt], '--', c='r')
ax1.set_xlim(minx, maxx)
ax1.set_ylim(0, 90)
ax1.set_xlabel(xlabel)
xt = ax1.get_xticks()
xt[xt < 0] += 24
ax1.set_xticklabels(['%i' % n for n in xt])
ax1.set_ylabel('Altitude [deg]')
ax1.set_title(titlenight, fontsize=11)

ax2 = ax1.twinx()
altitude = ax1.get_yticks()*u.deg
airmass = 1./np.cos(90*u.deg - altitude)
ax2.set_ylabel('Airmass')
myticks = []
for airval in airmass:
    if airval > 10:
        myticks.append('')
    else:
        myticks.append('%.2f' % airval)
ax2.set_yticklabels(myticks)
ax2.set_ylim(0, 90)
if not isList:
    plt.colorbar(sc, pad=.1).set_label('Azimuth [deg]')
plt.tight_layout()

if '--noshow' in sys.argv:
    if '--savechart' in sys.argv:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--savechart']
        try:
            givenname = np.array(sys.argv)[index + 1].item()
        except:
            raise IOError('%s is not a valid name for chart' % givenname)
        if (givenname == '') or (givenname.split('-')[0] == ''):
            raise IOError('blanck space is not a valid name for chart')
        else:
            if not os.path.isdir(os.getcwd() + '/figs'):
                os.mkdir(os.getcwd() + '/figs')
            figname = os.getcwd() + '/figs/' + givenname + '_' + night_beg + '_plan.png'
            print('map is being saved to', figname)
            plt.savefig(figname, format='png')
    plt.close()
else:
    if '--savechart' in sys.argv:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--savechart']
        try:
            givenname = np.array(sys.argv)[index + 1].item()
        except:
            raise IOError('blanck space is not a valid name for chart')
        if (givenname == '') or (givenname.split('-')[0] == ''):
            raise IOError('%s is not a valid name for chart' % givenname)
        else:
            if not os.path.isdir(os.getcwd() + '/figs'):
                os.mkdir(os.getcwd() + '/figs')
            figname = os.getcwd() + '/figs/' + givenname + '_' + night_beg + '_plan.png'
            print('map is being saved to', figname)
            plt.savefig(figname, format='png')
    plt.show(block=False)
