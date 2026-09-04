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


def track_summary(track, local_times, night_mask):
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

    Returns a dict with keys: id, name, swatch, ra, dec, blinit, blockdur,
    peakalt, peaktime, airmass, moondist -- all display strings except id,
    name and swatch.
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
        _hh, _mm = int(_dur_h), int(round((_dur_h - int(_dur_h)) * 60.))
        _blockdur_str = f'{_hh}h{_mm:02d}' if _hh else f'{_mm}m'
        _blinit_str = _hours_to_hhmm(track['block_starts'])
    else:
        _blockdur_str = '—'
        _blinit_str = '—'

    _moondist = track.get('moon_distance')
    _moondist_str = '%i' % _moondist if _moondist is not None else '—'

    return {'id': track['name'], 'name': track['name'],
           'swatch': track['color'], 'ra': _ra_str, 'dec': _dec_str,
           'blinit': _blinit_str, 'blockdur': _blockdur_str,
           'peakalt': _peakalt_str, 'peaktime': _peaktime_str,
           'airmass': _airmass_str, 'moondist': _moondist_str}


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
        if _t.get('has_block'):
            _blinit = _hours_to_hhmm(_t['block_starts']) + ':00'
            _blocktime = (_t['block_ends'] - _t['block_starts']) * 3600.
        else:
            _blinit = '00:00:00'
            _blocktime = 0.
        _lines.append('%s,%.6f,%.6f,%s,%.0f' % (
            _t['name'], _coords.ra.deg, _coords.dec.deg, _blinit,
            _blocktime))
    return '\n'.join(_lines) + '\n'
