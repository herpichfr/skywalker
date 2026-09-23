"""Pure helpers shared by the matplotlib and HTML renderers.

Nothing here touches matplotlib, plotly or astropy: it only reshapes arrays
and booleans that the Skywalker pipeline has already computed, so neither
renderer needs to duplicate any astronomy.
"""

from dataclasses import dataclass
from datetime import datetime

import numpy as np

# matplotlib's default colour cycle (tab10), so a track with no colour of its
# own (e.g. a track added outside set_plot()'s property-cycler loop) still
# matches the PNG's palette.
PALETTE = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

ASTRO_NIGHT_ALT_DEG = -18.0


def astro_night_mask(sun_alt):
    """True where the Sun is below -18 deg: astronomical night.

    This is a strictly tighter condition than the "Sun below the horizon"
    (0 deg) cut used by track_summary()'s night_mask argument and by
    cli.Skywalker.compute_track(): every sample where this mask is True is
    also True under that 0 deg cut, but the converse does not hold.
    """
    return np.asarray(sun_alt, dtype=float) < ASTRO_NIGHT_ALT_DEG


def palette_color(index):
    """Colour for the index-th object, cycling matplotlib's tab10 order."""
    return PALETTE[index % len(PALETTE)]


def airmass_from_alt(alt_deg):
    """Airmass at altitude(s) alt_deg, in degrees.

    Same formula as cli.Skywalker.set_plot() uses for the airmass twin axis:
    1/cos(90 - alt) = sec(z). Returns NaN below the horizon, where airmass is
    not meaningful.
    """
    alt = np.asarray(alt_deg, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        am = 1. / np.cos(np.deg2rad(90. - alt))
    return np.where(alt > 0., am, np.nan)


def mask_to_intervals(x, mask):
    """Return the [x_start, x_end] spans of contiguous True runs in mask.

    Parameters:
    -----------
    x : np.ndarray
        Coordinate array (e.g. decimal hours, or datetimes), same length as
        mask.
    mask : np.ndarray of bool
        Boolean array to find contiguous True runs in.
    """
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return []
    idx = np.flatnonzero(mask)
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([idx[0]], idx[breaks + 1]))
    ends = np.concatenate((idx[breaks], [idx[-1]]))
    return [(x[s], x[e]) for s, e in zip(starts, ends)]


@dataclass
class Band:
    """One twilight/darkness fill_between span of the altitude panel."""
    kind: str          # 'civil' | 'nautical' | 'astronomical' | 'dark' | 'moonlit'
    t0: datetime
    t1: datetime
    color: str
    alpha: float


def _hours_to_hhmm(hours):
    """Format decimal hours from local midnight as an HH:MM clock string.

    Same convention as hover.TimeCursor._clock_label(): negative hours are
    before midnight, so they wrap into the next day's clock time.
    """
    _h = (hours + 24. if hours < 0. else hours) % 24.
    _hh = int(_h)
    _mm = int(round((_h - _hh) * 60.))
    if _mm == 60:                       # e.g. 01:59:40 reads as 02:00
        _mm = 0
        _hh = (_hh + 1) % 24
    return '%02d:%02d' % (_hh, _mm)


def track_summary(track, local_times, night_mask, astro_mask=None,
                  minalt=None):
    """One DataTable row's worth of numbers for a track.

    Parameters:
    -----------
    track : dict
        A track registry entry, as built by cli.Skywalker.compute_track()
        (or the Moon's track, which is built separately and lacks 'coords').
    local_times : sequence of datetime.datetime
        Local wall-clock time of every sample, aligned with track['alt'].
    night_mask : np.ndarray of bool
        True where the Sun is below the horizon; the peak is searched only
        here, since a daytime culmination is useless for planning.
    astro_mask : np.ndarray of bool, optional
        True where the Sun is below -18 deg (astronomical night; see
        astro_night_mask()), aligned with track['alt']. Drives the
        HoursObs column, using the same night definition and the same
        (mask & alt > minalt).sum() * step_h computation as
        cli.Skywalker.year_max_altitudes()'s hours_up, so the two agree
        for the same date/target. None (the default) leaves HoursObs
        as '—', same as a track with no resolved coordinates.
    minalt : float, optional
        Minimum usable altitude in degrees, paired with astro_mask.

    Returns a dict with keys: id, name, swatch, ra, dec, blinit, blockdur,
    peakalt, peaktime, airmass, usablehours, moondist -- all display
    strings except id, name and swatch. blockdur is '%dh%02d' when the
    block is a whole number of hours and minutes and >= 1h (e.g. '1h30'),
    '%dm' when it is whole minutes under 1h (e.g. '45m'), and total
    seconds with an 's' suffix otherwise (e.g. '9015s'), so every value
    webapp._parse_blockdur() accepts reads back to the identical number
    of seconds -- that function inverts this formatting, and its
    docstring names the same forms. A track with no block gets '—'.
    """
    _alt = np.asarray(track['alt'], dtype=float)
    _night_alt = np.where(np.asarray(night_mask, dtype=bool), _alt, np.nan)
    if np.all(np.isnan(_night_alt)):
        _peakalt_str, _peaktime_str, _airmass_str = '—', '—', '—'
    else:
        _i = int(np.nanargmax(_night_alt))
        _peakalt = _night_alt[_i]
        _peakalt_str = '%.1f' % _peakalt
        _peaktime_str = local_times[_i].strftime('%H:%M')
        _air = airmass_from_alt(np.array([_peakalt]))[0]
        _airmass_str = '%.2f' % _air if np.isfinite(_air) else '—'

    _coords = track.get('coords')
    if _coords is not None:
        _ra_str = _coords.ra.to_string(unit='hourangle', sep=':',
                                       precision=1, pad=True)
        _dec_str = _coords.dec.to_string(unit='deg', sep=':', precision=0,
                                         alwayssign=True, pad=True)
    else:
        _ra_str = _dec_str = '—'

    if track.get('has_block'):
        _dur_h = track['block_ends'] - track['block_starts']
        # block_ends is derived from a seconds -> astropy-hour -> back-to-
        # seconds round trip, which leaves float noise (e.g.
        # 5399.999999999996) that would otherwise push a duration into the
        # wrong branch below or let _mm come out as 60. Rounding to the
        # nearest whole second recovers the exact integer every value
        # webapp._parse_blockdur() accepts is built from, which is what
        # makes this formatting round-trip through it exactly.
        _total_s = round(_dur_h * 3600.)
        _hh, _rem_s = divmod(_total_s, 3600)
        _mm, _ss = divmod(_rem_s, 60)
        if _ss == 0 and _hh >= 1:
            _blockdur_str = f'{_hh}h{_mm:02d}'
        elif _ss == 0:
            _blockdur_str = f'{_mm}m'
        else:
            _blockdur_str = f'{_total_s}s'
        _blinit_str = _hours_to_hhmm(track['block_starts'])
    else:
        _blockdur_str = '—'
        _blinit_str = '—'

    _moondist = track.get('moon_distance')
    _moondist_str = '%i' % _moondist if _moondist is not None else '—'

    if (astro_mask is not None and minalt is not None
            and _coords is not None):
        # Same definition and formula as
        # cli.Skywalker.year_max_altitudes()'s hours_up: astronomical
        # night (Sun below -18 deg) AND alt above minalt, summed and
        # scaled by the sample step in hours -- both grids are evenly
        # spaced over a 24h window, so step_h = 24 / (n_samples - 1)
        # holds for either one.
        _step_h = 24. / (len(_alt) - 1)
        _usable_mask = np.asarray(astro_mask, dtype=bool) & (_alt > minalt)
        _usablehours_str = '%.1f' % (_usable_mask.sum() * _step_h)
    else:
        _usablehours_str = '—'

    return {'id': track['name'], 'name': track['name'],
           'swatch': track['color'], 'ra': _ra_str, 'dec': _dec_str,
           'blinit': _blinit_str, 'blockdur': _blockdur_str,
           'peakalt': _peakalt_str, 'peaktime': _peaktime_str,
           'airmass': _airmass_str, 'usablehours': _usablehours_str,
           'moondist': _moondist_str}


def track_blinit_and_blocktime(track):
    """A track's block-start string and block length in seconds.

    Inverts the (block_starts, block_ends) hours pair that
    cli.Skywalker.compute_track() stores back into the "HH:MM:SS" +
    seconds form its own blinit/blocktime parameters take, so a track
    already held in memory can be fed straight back into compute_track()
    -- e.g. to recompute it at a new site, or after editing one of its
    other fields -- without re-deriving the user's original input. A
    track with no block returns ('00:00:00', 0.).
    """
    if not track.get('has_block'):
        return '00:00:00', 0.
    return (_hours_to_hhmm(track['block_starts']) + ':00',
           (track['block_ends'] - track['block_starts']) * 3600.)


def format_selection_csv(tracks):
    """Format tracks as a skywalker NAME,RA,DEC,BLINIT,BLOCKTIME CSV.

    RA/Dec are written as decimal degrees so the file round-trips under
    skywalker's default --raunit=auto: a bare decimal above 24 is always
    read as degrees (see coords.parse_ra()). The Moon and any track without
    resolved coordinates are silently skipped.
    """
    _lines = ['NAME,RA,DEC,BLINIT,BLOCKTIME']
    for _t in tracks:
        _coords = _t.get('coords')
        if _t.get('is_moon') or _coords is None:
            continue
        _blinit, _blocktime = track_blinit_and_blocktime(_t)
        _lines.append('%s,%.6f,%.6f,%s,%.0f' % (
            _t['name'], _coords.ra.deg, _coords.dec.deg, _blinit,
            _blocktime))
    return '\n'.join(_lines) + '\n'
