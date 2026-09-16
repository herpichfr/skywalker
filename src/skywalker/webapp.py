"""Interactive web UI: a target table below the plot, letting the user
include/remove targets and add arbitrary new ones (by name or coordinates)
while looking at the night. Uses Dash, an optional dependency -- importing
this module never requires it, only run_webapp() does.

NOTE: dash_table.DataTable is officially deprecated (removal planned for
Dash 5.0; the migration target is dash-ag-grid). pyproject.toml therefore
pins "dash>=2.14,<5" rather than an open-ended floor.

Architecture, to avoid a callback loop where several controls all want to
write sw-table.selected_rows while the figure reads it:
    - the MUTATOR callback (add/remove/all/none/invert) may only read
      data/selected_rows as State, and is the only writer of both.
    - the CONSUMER callback (redraws the figure) reads them as Input and
      writes nothing a mutator reads.
This makes the callback graph acyclic by construction, so no
allow_duplicate=True is needed anywhere.

State lives server-side in a TrackRegistry, never in dcc.Store: Dash
serialises Store values with plotly's JSON encoder, which silently turns
numpy arrays into lists and datetimes into strings rather than erroring, so
a track would come back as different Python types after a browser
round-trip. The registry is shared by every browser tab that connects to
this process (single-observer tool -- see --web's help text).
"""

import re
import threading

import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.coordinates.name_resolve import NameResolveError

from . import coords
from . import plotdata

_COLUMNS = [
    {'name': '', 'id': 'swatch'},
    {'name': 'Target', 'id': 'name', 'editable': True},
    {'name': 'RA', 'id': 'ra', 'editable': True},
    {'name': 'Dec', 'id': 'dec', 'editable': True},
    {'name': 'Block', 'id': 'blinit', 'editable': True},
    {'name': 'Dur', 'id': 'blockdur', 'editable': True},
    {'name': 'Peak alt', 'id': 'peakalt'},
    {'name': 'at', 'id': 'peaktime'},
    {'name': 'Best X', 'id': 'airmass'},
    {'name': 'Moon', 'id': 'moondist'},
]

_OK_STYLE = {'color': '#070', 'minHeight': '1.4em', 'fontSize': '13px'}
_ERR_STYLE = {'color': '#b00', 'minHeight': '1.4em', 'fontSize': '13px'}


def _require_dash():
    """Import dash lazily, raising a clear error if it is not installed."""
    try:
        from dash import (Dash, dcc, html, dash_table, Input, Output,
                          State, ctx, no_update)
    except ImportError as exc:
        raise ImportError(
            "--web requires dash, which is not installed.\n"
            "  Install it with:  pip install 'dash>=2.14,<5'\n"
            "  Or reinstall skywalker with the web extra:  "
            "pip install -e '.[web]'") from exc
    return Dash, dcc, html, dash_table, Input, Output, State, ctx, no_update


def _format_blinit(value):
    """Normalize a block-start time to HH:MM:SS.

    A local reimplementation rather than a call to
    cli.Skywalker.check_blockinit_format(): that method copies its input
    into a fixed-width numpy string array before padding, so a short value
    like "20:30" can silently truncate back to 5 characters after ":00" is
    appended. Plain string ops have no such width limit.

    Raises ValueError unless the result is a genuine 0 <= H < 24,
    0 <= M < 60, 0 <= S < 60 clock time, WITH one exemption: "24:00:00"
    (H==24 and M==S==0) is also accepted, because this same helper also
    normalizes a block END via _resolve_blocktime(), where "24:00" is
    standard end-of-day notation -- _hours_from_hhmmss() already
    normalizes "24:00:00" to 0.0, and _resolve_blocktime()'s past-midnight
    wrap then yields the correct duration from it, exactly as it does for
    "00:00:00". Refusing "24:00" here would silence that, not fix
    anything. This used to accept genuine garbage like "25:99:99"
    outright, silently, which is what the rest of this check still
    catches. All three call sites -- the Add-target form's block-start
    field, the table's editable Block cell, and the block-end field
    inside _resolve_blocktime() -- already wrap this call in a
    try/except, so raising here is enough to make all three reject bad
    input.
    """
    _value = str(value).strip()
    try:
        _dec = float(_value)
    except ValueError:
        _parts = _value.split(':')
        while len(_parts) < 3:
            _parts.append('00')
        _result = ':'.join(_parts)
    else:
        _h = int(_dec)
        _m = int(round((_dec - _h) * 60.))
        _result = f"{_h}:{_m:02d}:00"
    try:
        _h, _m, _s = (float(_p) for _p in _result.split(':'))
    except ValueError:
        raise ValueError(f"'{_value}' is not a valid HH:MM:SS time.")
    _is_end_of_day = (_h == 24. and _m == 0. and _s == 0.)
    if not _is_end_of_day and not (
            0. <= _h < 24. and 0. <= _m < 60. and 0. <= _s < 60.):
        raise ValueError(
            f"'{_value}' is not a valid HH:MM:SS time: hours must be "
            "0-23 (or exactly 24:00:00), minutes/seconds 0-59.")
    return _result


def _hours_from_hhmmss(hhmmss):
    """Decimal hours from local midnight for an "HH:MM:SS" string.

    Same wrap convention cli.Skywalker.compute_track() uses for
    block_starts, so a duration computed by subtracting two of these is
    consistent with block_starts/block_ends everywhere else.
    """
    _h, _m, _s = (float(_p) for _p in hhmmss.split(':'))
    _hours = _h + _m / 60. + _s / 3600.
    if _hours > 12.:
        _hours -= 24.
    return _hours


def _resolve_blocktime(blinit, blockend_text, blocktime_seconds):
    """Effective block length in seconds, or (None, error_message).

    blockend takes priority over blocktime_seconds when both are given (per
    user's explicit choice); with neither, the block length is 0, which
    cli.Skywalker.compute_track() treats as "no block" -- no shading under
    the altitude track.
    """
    _blockend_text = (blockend_text or '').strip()
    if _blockend_text:
        try:
            _blockend_fmt = _format_blinit(_blockend_text)
        except Exception as exc:
            return None, f"Bad block end '{_blockend_text}': {exc}"
        _duration_h = (_hours_from_hhmmss(_blockend_fmt)
                       - _hours_from_hhmmss(blinit))
        if _duration_h <= 0:
            return None, "Block ends must be after block starts."
        return _duration_h * 3600., None
    if blocktime_seconds:
        return float(blocktime_seconds), None
    return 0., None


def _sanitize_name(name):
    """Strip characters that would break a DataTable filter_query or the CSV."""
    return name.strip().replace('"', '').replace(',', '')


_BLOCKDUR_RE = re.compile(r'^(\d+)h(\d{1,2})$')


def _parse_blockdur(text):
    """Seconds from a Dur cell's "<H>h<MM>" / "<M>m" / "\u2014" display text.

    Inverts plotdata.track_summary()'s '%ih%02d'/'%im'/'\u2014' formatting,
    so editing the Dur cell round-trips through the same convention the
    table already displays. A bare number is read as seconds, matching
    the Add-target form's "Block size [s]" field. Raises ValueError on
    anything else.
    """
    _text = str(text).strip()
    if not _text or _text == '\u2014':
        return 0.
    try:
        return float(_text)
    except ValueError:
        pass
    _match = _BLOCKDUR_RE.match(_text)
    if _match:
        return int(_match.group(1)) * 3600. + int(_match.group(2)) * 60.
    if _text.endswith('m') and _text[:-1].isdigit():
        return float(_text[:-1]) * 60.
    raise ValueError(
        f"Bad block duration {_text!r}: use '1h30', '45m', or a number "
        "of seconds.")


class TrackRegistry:
    """Server-side store of the session's tracks, keyed by target name.

    Colours are assigned once, on first add, and never reassigned -- so
    removing and re-adding a target restores its original colour, and
    unchecking one target never recolours the others.
    """

    def __init__(self, walker):
        self.walker = walker
        self.tracks = {}            # name -> track dict, insertion-ordered
        self.colors = {}            # name -> hex colour, never pruned
        self.year_curves = {}       # name -> year curve dict
                                    # (see year_series()) -- a curve is a
                                    # pure function of a target's
                                    # coordinates and the site, so it goes
                                    # stale when either changes; a
                                    # coordinate edit pops just that
                                    # target's entry (see recompute()) and
                                    # a site change clears the whole cache
                                    # (see apply_site())
        self.color_cursor = 0
        self.lock = threading.Lock()
        self.local_times = (walker.frame_time_overnight.obstime
                            + walker.utcoffset).datetime
        self.night_mask = walker.sunaltaz_time_overnight.alt.value < 0.
        self.moon = {'name': 'Moon',
                     'alt': walker.moonaltaz_time_overnight.alt.value,
                     'az': walker.moonaltaz_time_overnight.az.value,
                     'color': 'c', 'is_moon': True, 'has_block': False,
                     'block_starts': 0., 'block_ends': 0.,
                     'moon_distance': None, 'label_x': None,
                     'label_y': None, 'label_color': None,
                     'chart_alt': None, 'chart_az': None,
                     'chart_hours': None, 'coords': None}

    def _next_color(self, name):
        if name not in self.colors:
            self.colors[name] = plotdata.palette_color(self.color_cursor)
            self.color_cursor += 1
        return self.colors[name]

    def add(self, name, ra_deg, dec_deg, blinit, blocktime=0.):
        """Compute and store one track. Returns (track, error_message)."""
        name = _sanitize_name(str(name))
        if not name:
            return None, "Enter a target name."
        with self.lock:
            if name in self.tracks:
                return None, f"'{name}' is already in the list."
            _color = self._next_color(name)
            _track = self.walker.compute_track(
                name, ra_deg, dec_deg, blinit, blocktime=blocktime,
                color=_color)
            if _track is None:
                return None, (
                    f"{name} never rises above the horizon while the Sun "
                    f"is down at {self.walker.sitename} on the night of "
                    f"{self.walker.nightstarts}.")
            self.tracks[name] = _track
            return _track, None

    def remove(self, names):
        """Remove tracks by name. The Moon can never be removed."""
        with self.lock:
            for _name in names:
                self.tracks.pop(_name, None)

    def recompute(self, name, ra_deg, dec_deg, blinit, blocktime=0.):
        """Recompute one existing target's track in place, same colour.

        Used by a cell edit (Feature 2) to turn a (possibly new)
        coordinate/block pair into a fresh track dict without
        duplicating cli.Skywalker.compute_track()'s call convention.
        NOT used by apply_site()'s per-target rebuild, which calls
        self.walker.compute_track() directly instead: this method takes
        self.lock, and apply_site() already holds it for the whole
        rebuild, so calling this from inside apply_site() would deadlock
        on self.lock, a non-reentrant threading.Lock.

        Returns (track, error_message). On error self.tracks is left
        untouched. Does not touch self.year_curves -- an edited
        coordinate goes stale, a re-timed block does not, so the caller
        decides whether to invalidate that target's curve.
        """
        with self.lock:
            if name not in self.tracks:
                return None, f"'{name}' is not a known target."
            _color = self.colors.get(name)
            _track = self.walker.compute_track(
                name, ra_deg, dec_deg, blinit, blocktime=blocktime,
                color=_color)
            if _track is None:
                return None, (
                    f"{name} never rises above the horizon while the Sun "
                    f"is down at {self.walker.sitename} on the night of "
                    f"{self.walker.nightstarts}.")
            self.tracks[name] = _track
            return _track, None

    def rename(self, old_name, new_name):
        """Rename an existing target, re-keying tracks/colors/year_curves.

        The rename keeps its colour and its cached year curve (if any) by
        moving the colors/year_curves entries to the new key instead of
        letting them be recomputed under it. Rejects a name already in
        use, including 'Moon', which can never be a target name.

        Returns (track, error_message). Leaves everything untouched on
        error.
        """
        _new_name = _sanitize_name(str(new_name))
        if not _new_name:
            return None, "Enter a target name."
        with self.lock:
            if old_name not in self.tracks:
                return None, f"'{old_name}' is not a known target."
            if _new_name == old_name:
                return self.tracks[old_name], None
            if _new_name in self.tracks or _new_name == 'Moon':
                return None, f"'{_new_name}' is already in the list."
            _new_tracks = {}
            for _n, _t in self.tracks.items():
                if _n == old_name:
                    _t['name'] = _new_name
                    _new_tracks[_new_name] = _t
                else:
                    _new_tracks[_n] = _t
            self.tracks = _new_tracks
            if old_name in self.colors:
                self.colors[_new_name] = self.colors.pop(old_name)
            if old_name in self.year_curves:
                self.year_curves[_new_name] = self.year_curves.pop(old_name)
            return self.tracks[_new_name], None

    def ordered_tracks(self, names=None):
        """Tracks in insertion order (Moon last), optionally filtered."""
        _tracks = list(self.tracks.values())
        if names is not None:
            _keep = set(names)
            _tracks = [t for t in _tracks if t['name'] in _keep]
        _include_moon = names is None or 'Moon' in set(names)
        return _tracks + ([self.moon] if _include_moon else [])

    def table_rows(self):
        """One plotdata.track_summary() row per track, Moon included."""
        return [plotdata.track_summary(t, self.local_times, self.night_mask)
               for t in self.ordered_tracks()]

    def year_series(self, names=None):
        """Cached per-target year curves, computed on first use.

        Keyed by name and pruned only when it must be: a target's curve
        is a pure function of its coordinates and the site, so removing
        and re-adding a target (remove() deliberately has no matching
        pop() for this cache) gets its curve back for free, but an
        edited ra/dec (see recompute() callers) pops just that target's
        entry, and a site change (see apply_site()) clears the whole
        cache -- callers elsewhere must not assume an entry is
        permanent.

        The Moon, and any track with no resolved coordinates, is skipped: a
        monthly sample of the Moon's peak altitude is synodic aliasing, not
        information.
        """
        self.walker.set_year_frames()
        _tracks = [t for t in self.ordered_tracks(names)
                  if not t.get('is_moon') and t.get('coords') is not None]
        _missing = [t for t in _tracks if t['name'] not in self.year_curves]
        if _missing:
            _ra_deg = np.array([t['coords'].ra.deg for t in _missing])
            _dec_deg = np.array([t['coords'].dec.deg for t in _missing])
            _obj_coords = SkyCoord(ra=_ra_deg, dec=_dec_deg,
                                   unit=('deg', 'deg'))
            _results = self.walker.year_max_altitudes(_obj_coords)
            with self.lock:
                for _t, _res in zip(_missing, _results):
                    self.year_curves[_t['name']] = _res
        _curves = [dict(self.year_curves[t['name']], name=t['name'],
                       color=self.colors.get(t['name']))
                  for t in _tracks]
        return _curves, self.walker.year_dates

    def apply_site(self, site=None, lat=None, lon=None, elev=None,
                   name=None):
        """Move the whole session to a new observing site, in place.

        lat and lon, given together, take precedence over site -- this
        mirrors the Lat/Lon inputs overriding the site dropdown in the
        web UI. Rebuilds the walker's location, observer, time frames and
        night frames by calling cli.Skywalker's own set_location/
        set_observer/set_time/set_night_frames (never duplicating their
        logic), invalidates the year-view cache, and recomputes every
        stored target's track at the new site from its saved coordinates.
        Colours are preserved: self.colors is never touched here.

        name, given only alongside lat/lon, becomes walker.sitename in
        place of the "{lat} {lon}" fallback -- see
        cli.Skywalker.set_location()'s own name parameter, which this
        passes straight through and which does the actual stripping and
        fallback. Ignored when site is given instead.

        Returns (dropped_names, error_message). dropped_names lists any
        target that no longer rises above the horizon at the new site --
        it is removed from self.tracks, never left in a half-computed
        state. On error, nothing is left changed: the lat/lon range check
        runs before any attribute is mutated, and every walker attribute
        touched afterwards is snapshotted first and restored if a later
        step raises -- e.g. set_time()'s timezone lookup failing over
        open ocean.
        """
        _walker = self.walker
        with self.lock:
            _use_latlon = lat is not None and lon is not None
            if _use_latlon and not (-180. <= lon <= 360.):
                return [], f"Longitude {lon} out of range (-180 to 360 deg)."

            _snapshot = {_attr: getattr(_walker, _attr) for _attr in (
                'location', 'sitename', 'observer', 'inithour', 'utcoffset',
                'obs_time', 'delta_midnight', 'frame_time_overnight',
                'moon', 'sunaltaz_time_overnight', 'moonaltaz_time_overnight',
                'moon_brightness', 'year_frame', 'year_dates', 'year_shape',
                'year_night_mask', 'year_local_times')}
            try:
                if _use_latlon:
                    _walker.set_location(lat=lat, lon=lon, elev=elev,
                                        name=name)
                else:
                    _walker.set_location(site=site)
                _walker.set_observer()
                _walker.set_time()
                _walker.set_night_frames()
            except Exception as exc:
                for _attr, _val in _snapshot.items():
                    setattr(_walker, _attr, _val)
                return [], f"Could not apply site: {exc}"

            # Committed past this point: force set_year_frames() to
            # rebuild on next use, and drop every cached year curve --
            # both are pure functions of coordinates AND site now.
            _walker.year_frame = None
            _walker.year_dates = None
            _walker.year_shape = None
            _walker.year_night_mask = None
            _walker.year_local_times = None
            self.year_curves = {}

            self.local_times = (_walker.frame_time_overnight.obstime
                                + _walker.utcoffset).datetime
            self.night_mask = _walker.sunaltaz_time_overnight.alt.value < 0.
            self.moon = {'name': 'Moon',
                        'alt': _walker.moonaltaz_time_overnight.alt.value,
                        'az': _walker.moonaltaz_time_overnight.az.value,
                        'color': 'c', 'is_moon': True, 'has_block': False,
                        'block_starts': 0., 'block_ends': 0.,
                        'moon_distance': None, 'label_x': None,
                        'label_y': None, 'label_color': None,
                        'chart_alt': None, 'chart_az': None,
                        'chart_hours': None, 'coords': None}

            _dropped = []
            _new_tracks = {}
            for _name, _track in self.tracks.items():
                _coords = _track.get('coords')
                if _coords is None:
                    continue
                _blinit, _blocktime = plotdata.track_blinit_and_blocktime(
                    _track)
                _recomputed = _walker.compute_track(
                    _name, _coords.ra.deg, _coords.dec.deg, _blinit,
                    blocktime=_blocktime, color=self.colors.get(_name))
                if _recomputed is None:
                    _dropped.append(_name)
                    continue
                _new_tracks[_name] = _recomputed
            self.tracks = _new_tracks
            return _dropped, None


def _resolve_target(registry, name, ra_text, dec_text, raunit):
    """Resolve a name/RA/Dec triple to (ra_deg, dec_deg, name, error).

    RA+Dec, when both given, are read locally with no network. Otherwise
    the name is looked up via Sesame, which needs a network connection --
    the failure modes are distinguished so an offline user gets an
    actionable message instead of a bare exception.
    """
    _name = _sanitize_name(str(name or '').strip())
    _ra_text = (ra_text or '').strip()
    _dec_text = (dec_text or '').strip()

    if _ra_text and _dec_text:
        try:
            _ra = coords.parse_ra(_ra_text, raunit=raunit)
            _dec = coords.parse_dec(_dec_text)
        except ValueError as exc:
            return None, None, None, str(exc)
        if not _name:
            _name = f"Obj{len(registry.tracks) + 1}"
        return _ra, _dec, _name, None

    if _ra_text or _dec_text:
        return None, None, None, "Give both RA and Dec (or just a name)."

    if not _name:
        return None, None, None, "Enter an object name, or an RA and a Dec."

    try:
        _coord = SkyCoord.from_name(_name)
    except NameResolveError as exc:
        if 'All Sesame queries failed' in str(exc):
            return None, None, None, (
                "No network: Sesame name lookup failed. "
                "Enter RA and Dec directly.")
        return None, None, None, (
            f"Could not resolve '{_name}'. Check the spelling, "
            "or enter RA and Dec directly.")
    except Exception as exc:                       # pragma: no cover
        return None, None, None, f"Name lookup failed: {exc}"
    return _coord.ra.deg, _coord.dec.deg, _name, None


def _swatch_styles(rows, focus=None):
    """style_data_conditional: one colour chip per distinct palette colour,
    plus a highlight for the focused row. One condition per colour, not per
    row, so this stays small (<=10) no matter how many targets are listed.
    """
    _styles = []
    for _hex in sorted({r['swatch'] for r in rows if r.get('swatch')}):
        _styles.append({'if': {'filter_query': f'{{swatch}} = "{_hex}"',
                               'column_id': 'swatch'},
                        'backgroundColor': _hex, 'color': _hex})
    if focus:
        _styles.append({'if': {'filter_query': f'{{name}} = "{focus}"'},
                        'backgroundColor': '#eef4ff', 'fontWeight': 'bold'})
    return _styles


def _apply_focus(fig, focus):
    """Dim every trace/annotation that does not belong to the focused target."""
    if not focus:
        return fig
    for _trace in fig.data:
        _mine = (_trace.name == focus)
        _trace.opacity = 1.0 if _mine else 0.15
        if getattr(_trace, 'mode', None) == 'lines':
            _trace.line.width = 3 if _mine else 1
    for _ann in fig.layout.annotations or ():
        _ann.opacity = 1.0 if _ann.name == focus else 0.15
    return fig


def _resolve_focus(active_cell, rows):
    """active_cell's row_id, or None if it names no row currently in rows.

    _view and _view_year both read active_cell.get('row_id') as the name
    to focus on, but active_cell is client-side state that a rename, a
    removal, or a target dropped by a site switch never updates or
    clears -- it can go on naming a target that no longer exists. Left
    unchecked, that stale name matches no trace in _apply_focus() (or no
    row in _swatch_styles()), so every trace/annotation dims and the
    whole figure looks broken. Treating a focus that names no current
    row as "no focus" degrades to the correct, unfocused view instead of
    a duplicate Output on sw-table.active_cell (which would collide with
    _clear_highlight's) or a second writer next to the merged _mutate
    callback.
    """
    if not active_cell:
        return None
    _focus = active_cell.get('row_id')
    _names = {r['name'] for r in rows}
    return _focus if _focus in _names else None


_INPUT_STYLE = {'fontSize': '14px', 'padding': '5px 7px',
                'height': '30px', 'boxSizing': 'border-box'}


def _input_style(width):
    return dict(_INPUT_STYLE, width=width)


def _field(dcc, html, label, component):
    return html.Div([html.Label(label, style={'display': 'block',
                                              'fontSize': '11px',
                                              'color': '#555'}),
                     component],
                    style={'display': 'flex', 'flexDirection': 'column'})


def build_app(walker):
    """Build the Dash app for walker's already-computed night.

    Parameters:
    -----------
    walker : cli.Skywalker
        Must already have run set_location/set_observer/set_time/
        set_night_frames/set_target_list (main() does this before calling
        set_webapp()).
    """
    (Dash, dcc, html, dash_table, Input, Output, State, ctx,
     no_update) = _require_dash()
    from . import htmlplot

    registry = TrackRegistry(walker)
    for _idx in walker.target_list.index:
        _row = walker.target_list.loc[_idx]
        _, _err = registry.add(_row['NAME'], _row['RA'], _row['DEC'],
                               _row['BLINIT'], blocktime=_row['BLOCKTIME'])
        if _err is not None:
            walker.logger.warning(f"Skipping {_row['NAME']}: {_err}")

    # Computed once at startup: EarthLocation.get_site_names() hits
    # astropy's (possibly network-backed) site registry, so the dropdown's
    # option list is not refetched on every browser page load or callback.
    _site_options = [{'label': _s, 'value': _s}
                     for _s in sorted(set(EarthLocation.get_site_names()))
                     if _s]

    def _build_figure(names, focus=None):
        _fig = htmlplot.build_figure(
            registry.ordered_tracks(names), registry.local_times,
            walker.delta_midnight.value, walker.sunaltaz_time_overnight.alt.value,
            walker.moonaltaz_time_overnight.alt.value,
            walker.moon_brightness.value, walker.minalt,
            int(walker.utcoffset.value), walker.sitename, walker.nightstarts,
            make_skychart=walker.make_skychart)
        _fig.update_layout(width=None, autosize=True)
        return _apply_focus(_fig, focus)

    def _serve_layout():
        _rows = registry.table_rows()
        return html.Div([
            html.H3(f"Night starts: {walker.nightstarts} @ {walker.sitename}",
                   id='sw-title'),
            dcc.Loading(children=[
                dcc.Graph(id='sw-graph', figure=_build_figure(None),
                          config={'displaylogo': False,
                                  'modeBarButtonsToRemove': ['lasso2d',
                                                             'select2d']},
                          style={'width': '100%', 'height': '620px'}),
            ], type='default', delay_show=300),
            html.Button('Show year view', id='sw-year-btn', n_clicks=0),
            # Bare boolean -- safe in a dcc.Store. The module docstring's
            # warning is about numpy arrays and datetimes silently mangled
            # by plotly's JSON encoder; a plain bool round-trips exactly.
            dcc.Store(id='sw-year-on', data=False),
            html.Div([
                dcc.RadioItems(
                    id='sw-year-metric',
                    options=[{'label': 'Peak altitude', 'value': 'alt'},
                            {'label': 'Hours usable', 'value': 'hours'}],
                    value='alt', inline=True,
                    style={'fontSize': '13px', 'margin': '4px 0'}),
                dcc.Loading(children=[
                    dcc.Graph(id='sw-year', figure={},
                             config={'displaylogo': False},
                             style={'width': '100%', 'height': '380px'}),
                ], type='default', delay_show=300),
            ], id='sw-year-wrap', style={'display': 'none'}),
            html.Div([
                html.Span('Observatory', style={'fontWeight': 'bold',
                                                'fontSize': '13px',
                                                'alignSelf': 'center'}),
                _field(dcc, html, 'Site', dcc.Dropdown(
                    id='sw-site-dropdown', options=_site_options,
                    value=None, searchable=True, clearable=True,
                    placeholder='Search known sites...',
                    style={'width': '220px', 'fontSize': '13px'})),
                _field(dcc, html, 'Lat', dcc.Input(
                    id='sw-site-lat', type='number', placeholder='deg',
                    style=_input_style('90px'))),
                _field(dcc, html, 'Lon', dcc.Input(
                    id='sw-site-lon', type='number', placeholder='deg',
                    style=_input_style('90px'))),
                _field(dcc, html, 'Elev', dcc.Input(
                    id='sw-site-elev', type='number', placeholder='m',
                    style=_input_style('90px'))),
                _field(dcc, html, 'Name (optional)', dcc.Input(
                    id='sw-site-name', type='text',
                    placeholder='Display name', style=_input_style('160px'))),
                html.Button('Apply site', id='sw-site-apply', n_clicks=0),
            ], id='sw-site-row', style={'display': 'flex', 'gap': '8px',
                                        'alignItems': 'flex-end',
                                        'flexWrap': 'wrap',
                                        'margin': '10px 0'}),
            html.Div([
                _field(dcc, html, 'Name', dcc.Input(
                    id='sw-in-name', type='text', debounce=True,
                    placeholder='M31', style=_input_style('120px'))),
                _field(dcc, html, 'RA', dcc.Input(
                    id='sw-in-ra', type='text', debounce=True,
                    placeholder='00:42:44 or 10.68',
                    style=_input_style('170px'))),
                _field(dcc, html, 'Dec', dcc.Input(
                    id='sw-in-dec', type='text', debounce=True,
                    placeholder='+41:16:09 or 41.27',
                    style=_input_style('170px'))),
                _field(dcc, html, 'Block starts', dcc.Input(
                    id='sw-in-blinit', type='text',
                    placeholder=walker.time, style=_input_style('110px'))),
                _field(dcc, html, 'Block ends (optional)', dcc.Input(
                    id='sw-in-blockend', type='text', placeholder='—',
                    style=_input_style('110px'))),
                _field(dcc, html, 'Block size [s] (optional)', dcc.Input(
                    id='sw-in-blocktime', type='number', min=0, step=60,
                    style=_input_style('130px'))),
                html.Button('Add target', id='sw-add', n_clicks=0),
                html.Button('Remove selected', id='sw-remove', n_clicks=0),
                html.Span('|'),
                html.Button('All', id='sw-all', n_clicks=0),
                html.Button('None', id='sw-none', n_clicks=0),
                html.Button('Invert', id='sw-invert', n_clicks=0),
                html.Span('|'),
                html.Button('Clear highlight', id='sw-unhighlight',
                           n_clicks=0),
                dcc.Clipboard(id='sw-clip',
                             title='Copy selection as skywalker CSV '
                                   '(decimal-degree RA/Dec)',
                             style={'fontSize': '20px', 'cursor': 'pointer'}),
                html.Button('Download CSV', id='sw-dl-btn', n_clicks=0),
            ], id='sw-controls', style={'display': 'flex', 'gap': '8px',
                                        'alignItems': 'flex-end',
                                        'flexWrap': 'wrap',
                                        'margin': '10px 0'}),
            html.Div(id='sw-status', style=_OK_STYLE),
            dash_table.DataTable(
                id='sw-table', columns=_COLUMNS, data=_rows,
                selected_rows=list(range(len(_rows))),
                row_selectable='multi', sort_action='native',
                cell_selectable=True, page_action='none',
                style_table={'width': '100%', 'overflowX': 'auto',
                            'maxHeight': '340px', 'overflowY': 'auto'},
                style_cell={'fontFamily': 'monospace', 'fontSize': '12px',
                           'padding': '4px 8px', 'textAlign': 'right'},
                style_cell_conditional=[
                    {'if': {'column_id': 'swatch'}, 'width': '26px',
                     'maxWidth': '26px', 'overflow': 'hidden',
                     'padding': '4px 0'},
                    {'if': {'column_id': 'name'}, 'textAlign': 'left'}],
                style_header={'fontWeight': 'bold',
                             'backgroundColor': '#f2f2f2'},
                style_data_conditional=_swatch_styles(_rows),
                fixed_rows={'headers': True}),
            dcc.Download(id='sw-dl'),
        ], style={'maxWidth': '1420px', 'margin': '0 auto',
                 'fontFamily': 'sans-serif'})

    app = Dash(__name__, title=f"skywalker {walker.nightstarts}",
              update_title=None)
    app.layout = _serve_layout

    # Both the site switch (Feature 1) and per-cell edits (Feature 2) need
    # to write sw-table.data/selected_rows and sw-status -- the same
    # Outputs the MUTATOR above already owns -- so both are folded into
    # this one callback as two more ctx.triggered_id branches, rather than
    # adding a second writer and needing allow_duplicate=True. sw-title is
    # a new Output added to the same callback for the same reason: only
    # the site-apply branch has anything to say about it, every other
    # branch returns no_update for it.
    #
    # Editable-column edits arrive as Input('sw-table', 'data_timestamp'),
    # which dash_table only bumps for a front-end (user) edit of the
    # data prop -- never for a plain server-side Output write to data --
    # so this callback's own writes to sw-table.data cannot re-trigger the
    # 'sw-table' branch below and loop.
    _EDIT_FIELDS = ('name', 'ra', 'dec', 'blinit', 'blockdur')

    def _apply_cell_edit(old_row, field, new_value, selected):
        _name = old_row['name']
        if _name == 'Moon':
            return (registry.table_rows(), selected,
                   "The Moon row cannot be edited.", _ERR_STYLE, no_update,
                   no_update)
        _track = None
        try:
            if field == 'name':
                _track, _err = registry.rename(_name, new_value)
            elif field in ('ra', 'dec'):
                # Parse only the field that actually changed; take the
                # other one from the stored track in degrees. Re-parsing
                # BOTH via _resolve_target would feed it old_row['ra']/
                # ['dec'] -- plotdata.track_summary()'s display-formatted,
                # unit-less sexagesimal strings -- back through
                # coords.parse_ra(), which correctly refuses an
                # un-suffixed sexagesimal RA as ambiguous unless raunit
                # was explicitly set. That is the bug this avoids: the
                # unedited field never needs parsing at all.
                _old_track = registry.tracks[_name]
                try:
                    if field == 'ra':
                        _ra_deg = coords.parse_ra(
                            new_value, raunit=walker.raunit,
                            context=f" for {_name}")
                        _dec_deg = _old_track['coords'].dec.deg
                    else:
                        _ra_deg = _old_track['coords'].ra.deg
                        _dec_deg = coords.parse_dec(new_value)
                except ValueError as exc:
                    _err = str(exc)
                else:
                    _blinit, _blocktime = plotdata.track_blinit_and_blocktime(
                        _old_track)
                    _track, _err = registry.recompute(
                        _name, _ra_deg, _dec_deg, _blinit,
                        blocktime=_blocktime)
                    if _err is None:
                        with registry.lock:
                            registry.year_curves.pop(_name, None)
            elif field == 'blinit':
                _old_track = registry.tracks[_name]
                _, _blocktime = plotdata.track_blinit_and_blocktime(
                    _old_track)
                try:
                    _blinit = _format_blinit(new_value)
                except Exception as exc:
                    _err = f"Bad block start '{new_value}': {exc}"
                else:
                    _track, _err = registry.recompute(
                        _name, _old_track['coords'].ra.deg,
                        _old_track['coords'].dec.deg, _blinit,
                        blocktime=_blocktime)
            else:                                    # field == 'blockdur'
                _old_track = registry.tracks[_name]
                _blinit, _ = plotdata.track_blinit_and_blocktime(_old_track)
                try:
                    _seconds = _parse_blockdur(new_value)
                except ValueError as exc:
                    _err = str(exc)
                else:
                    _blocktime, _err = _resolve_blocktime(_blinit, None,
                                                          _seconds)
                    if _err is None:
                        _track, _err = registry.recompute(
                            _name, _old_track['coords'].ra.deg,
                            _old_track['coords'].dec.deg, _blinit,
                            blocktime=_blocktime)
        except Exception as exc:            # pragma: no cover -- safety net
            _track, _err = None, f"Bad edit to {_name}: {exc}"

        if _err is not None:
            return (registry.table_rows(), selected, _err, _ERR_STYLE,
                   no_update, no_update)
        return (registry.table_rows(), selected, f"Updated {_track['name']}.",
               _OK_STYLE, no_update, no_update)

    def _handle_cell_edit(rows, selected):
        _canon = registry.table_rows()
        if len(rows) != len(_canon):
            # A row appeared/disappeared without going through the
            # mutator above -- should not happen; resync defensively.
            return (_canon, list(range(len(_canon))), '', _OK_STYLE,
                   no_update, no_update)
        for _old, _new in zip(_canon, rows):
            for _field in _EDIT_FIELDS:
                if _old.get(_field) != _new.get(_field):
                    return _apply_cell_edit(_old, _field, _new[_field],
                                            selected)
        # No field actually differs (e.g. a cell was focused and blurred
        # unchanged): '' is safe here, not no_update, because dash_table's
        # installed async-table.js only stamps data_timestamp inside its
        # own __setProps wrapper -- i.e. when the table's own front-end
        # edit handling calls its outgoing setProps({data: ...}) to push
        # a user edit up to Dash. A callback Output writing sw-table.data
        # updates the prop top-down via React and never goes through that
        # wrapper, so this callback's own writes cannot re-trigger this
        # 'sw-table' branch and there is no risk of this '' clobbering a
        # freshly-written "Updated X." message.
        return rows, selected, '', _OK_STYLE, no_update, no_update

    @app.callback(
        Output('sw-table', 'data'),
        Output('sw-table', 'selected_rows'),
        Output('sw-status', 'children'),
        Output('sw-status', 'style'),
        Output('sw-in-name', 'value'),
        Output('sw-title', 'children'),
        Input('sw-add', 'n_clicks'),
        Input('sw-remove', 'n_clicks'),
        Input('sw-all', 'n_clicks'),
        Input('sw-none', 'n_clicks'),
        Input('sw-invert', 'n_clicks'),
        Input('sw-in-name', 'n_submit'),
        Input('sw-in-dec', 'n_submit'),
        Input('sw-site-apply', 'n_clicks'),
        Input('sw-table', 'data_timestamp'),
        State('sw-table', 'data'),
        State('sw-table', 'selected_rows'),
        State('sw-in-name', 'value'),
        State('sw-in-ra', 'value'),
        State('sw-in-dec', 'value'),
        State('sw-in-blinit', 'value'),
        State('sw-in-blockend', 'value'),
        State('sw-in-blocktime', 'value'),
        State('sw-site-dropdown', 'value'),
        State('sw-site-lat', 'value'),
        State('sw-site-lon', 'value'),
        State('sw-site-elev', 'value'),
        State('sw-site-name', 'value'),
        prevent_initial_call=True)
    def _mutate(_add, _remove, _all, _none, _invert, _sub1, _sub2, _apply,
               _ts, rows, selected, name, ra, dec, blinit, blockend,
               blocktime, site, lat, lon, elev, site_name):
        who = ctx.triggered_id
        rows = rows or []
        selected = selected or []

        if who in ('sw-add', 'sw-in-name', 'sw-in-dec'):
            _ra_deg, _dec_deg, _name, _err = _resolve_target(
                registry, name, ra, dec, walker.raunit)
            if _err is not None:
                return rows, selected, _err, _ERR_STYLE, name, no_update
            _blinit = blinit.strip() if blinit and blinit.strip() \
                else walker.time
            try:
                _blinit = _format_blinit(_blinit)
            except Exception as exc:
                return (rows, selected,
                       f"Bad block start '{_blinit}': {exc}", _ERR_STYLE,
                       name, no_update)
            _blocktime, _err = _resolve_blocktime(_blinit, blockend,
                                                  blocktime)
            if _err is not None:
                return rows, selected, _err, _ERR_STYLE, name, no_update
            _track, _err = registry.add(_name, _ra_deg, _dec_deg, _blinit,
                                        blocktime=_blocktime)
            if _err is not None:
                return rows, selected, _err, _ERR_STYLE, name, no_update
            _rows = registry.table_rows()
            _selected = list(range(len(_rows)))   # newly-added stays checked
            return (_rows, _selected,
                   f"Added {_track['name']} ({len(registry.tracks)} "
                   "target(s)).", _OK_STYLE, '', no_update)

        if who == 'sw-remove':
            _names = [rows[i]['name'] for i in selected
                     if i < len(rows) and rows[i]['name'] != 'Moon']
            registry.remove(_names)
            _rows = registry.table_rows()
            return (_rows, list(range(len(_rows))),
                   f"Removed {len(_names)} target(s).", _OK_STYLE, no_update,
                   no_update)

        if who == 'sw-all':
            return (rows, list(range(len(rows))), '', _OK_STYLE, no_update,
                   no_update)
        if who == 'sw-none':
            return rows, [], '', _OK_STYLE, no_update, no_update
        if who == 'sw-invert':
            _sel = set(selected)
            return (rows, [i for i in range(len(rows)) if i not in _sel],
                   '', _OK_STYLE, no_update, no_update)

        if who == 'sw-site-apply':
            _dropped, _err = registry.apply_site(
                site=site or None,
                lat=lat if lat not in (None, '') else None,
                lon=lon if lon not in (None, '') else None,
                elev=elev if elev not in (None, '') else None,
                name=site_name)
            if _err is not None:
                return rows, selected, _err, _ERR_STYLE, no_update, no_update
            _rows = registry.table_rows()
            _msg = f"Site set to {walker.sitename}."
            if _dropped:
                _msg += (" No longer observable, removed: "
                        + ', '.join(_dropped) + '.')
            _title = (f"Night starts: {walker.nightstarts} @ "
                     f"{walker.sitename}")
            return (_rows, list(range(len(_rows))), _msg, _OK_STYLE,
                   no_update, _title)

        if who == 'sw-table':
            return _handle_cell_edit(rows, selected)

        return rows, selected, '', _OK_STYLE, no_update, no_update

    # Exclusivity between the two site-input modes: selecting a dropdown
    # site clears Lat/Lon/Elev/Name, and typing a Lat/Lon/Elev value
    # clears the dropdown. A single callback with every site control as
    # both Input and Output -- Dash permits an Output that is also its
    # own Input; it is a cycle SPLIT ACROSS two callbacks that Dash
    # refuses at registration. The guard on the triggering control's own
    # new value being empty is what makes this converge instead of
    # oscillating: a dropdown pick clears lat/lon/elev, each of which
    # refires this same callback with its own (now empty) value as the
    # trigger, and the guard turns that refire into a no-op rather than a
    # bounce that wipes the dropdown pick that just landed. Name is an
    # Output only, never an Input, here -- typing a name must not clear
    # the dropdown.
    @app.callback(
        Output('sw-site-dropdown', 'value'),
        Output('sw-site-lat', 'value'),
        Output('sw-site-lon', 'value'),
        Output('sw-site-elev', 'value'),
        Output('sw-site-name', 'value'),
        Input('sw-site-dropdown', 'value'),
        Input('sw-site-lat', 'value'),
        Input('sw-site-lon', 'value'),
        Input('sw-site-elev', 'value'),
        prevent_initial_call=True)
    def _enforce_site_exclusivity(site, lat, lon, elev):
        who = ctx.triggered_id
        _triggered_value = {'sw-site-dropdown': site, 'sw-site-lat': lat,
                            'sw-site-lon': lon,
                            'sw-site-elev': elev}[who]
        if _triggered_value in (None, ''):
            return no_update, no_update, no_update, no_update, no_update
        if who == 'sw-site-dropdown':
            return no_update, None, None, None, ''
        return None, no_update, no_update, no_update, no_update

    @app.callback(
        Output('sw-graph', 'figure'),
        Output('sw-table', 'style_data_conditional'),
        Output('sw-clip', 'content'),
        Input('sw-table', 'data'),
        Input('sw-table', 'selected_rows'),
        Input('sw-table', 'active_cell'))
    def _view(rows, selected, active_cell):
        rows = rows or []
        selected = selected or []
        _names = [rows[i]['name'] for i in selected if i < len(rows)]
        _focus = _resolve_focus(active_cell, rows)
        _fig = _build_figure(_names, focus=_focus)
        _csv = plotdata.format_selection_csv(
            registry.ordered_tracks(_names))
        return _fig, _swatch_styles(rows, focus=_focus), _csv

    @app.callback(
        Output('sw-year-on', 'data'),
        Output('sw-year-wrap', 'style'),
        Output('sw-year-btn', 'children'),
        Input('sw-year-btn', 'n_clicks'),
        State('sw-year-on', 'data'),
        prevent_initial_call=True)
    def _toggle_year(_n, on):
        on = not on
        _style = {'display': 'block'} if on else {'display': 'none'}
        _label = 'Hide year view' if on else 'Show year view'
        return on, _style, _label

    @app.callback(
        Output('sw-year', 'figure'),
        Input('sw-year-on', 'data'),
        Input('sw-year-metric', 'value'),
        Input('sw-table', 'data'),
        Input('sw-table', 'selected_rows'),
        Input('sw-table', 'active_cell'))
    def _view_year(on, metric, rows, selected, active_cell):
        # Lazy contract: while the panel is hidden nothing here is
        # computed, so app startup and ordinary target-adding stay exactly
        # as fast as they are today.
        if not on:
            return no_update
        rows = rows or []
        selected = selected or []
        _names = [rows[i]['name'] for i in selected if i < len(rows)]
        _focus = _resolve_focus(active_cell, rows)
        _curves, _dates = registry.year_series(_names)
        if not _curves:
            return {}
        _fig = htmlplot.build_year_figure(
            _curves, _dates, walker.minalt, walker.sitename,
            walker.nightstarts, metric=metric)
        _fig.update_layout(width=None, autosize=True)
        return _apply_focus(_fig, _focus)

    @app.callback(Output('sw-table', 'active_cell'),
                 Input('sw-unhighlight', 'n_clicks'),
                 prevent_initial_call=True)
    def _clear_highlight(_n):
        return None

    @app.callback(Output('sw-dl', 'data'),
                 Input('sw-dl-btn', 'n_clicks'),
                 State('sw-table', 'data'), State('sw-table', 'selected_rows'),
                 prevent_initial_call=True)
    def _download(_n, rows, selected):
        rows = rows or []
        selected = selected or []
        _names = [rows[i]['name'] for i in selected if i < len(rows)]
        _csv = plotdata.format_selection_csv(registry.ordered_tracks(_names))
        return dcc.send_string(_csv, f"skywalker_{walker.nightstarts}.csv")

    return app


_LOOPBACK = {'127.0.0.1', 'localhost', '::1'}


def run_webapp(walker, host='127.0.0.1', port=8050, debug=False,
              open_browser=False):
    """Build and serve the interactive web app for walker's night.

    Parameters:
    -----------
    walker : cli.Skywalker
        Must already have completed set_location through set_target_list.
    host : str, optional
        Interface to bind. Default '127.0.0.1' (loopback only). A
        non-loopback host has no authentication in front of it -- use an
        SSH tunnel (ssh -L 8050:localhost:8050 user@host) for remote access
        instead, unless you understand and accept the exposure.
    port : int, optional
        TCP port. Default 8050.
    debug : bool, optional
        Enable Dash/Werkzeug's debugger. Refused when host is not
        loopback, since the debugger allows remote code execution.
    open_browser : bool, optional
        Open the UI in the default browser once the server is listening.
    """
    if debug and host not in _LOOPBACK:
        raise ValueError(
            f"Refusing to combine --web-debug with --web-host {host}: "
            "the interactive debugger allows remote code execution. Use "
            "a loopback host, or drop --web-debug.")
    if host not in _LOOPBACK:
        walker.logger.warning(
            f"Binding to {host} exposes an unauthenticated app to the "
            "network; prefer ssh -L 8050:localhost:8050 user@host.")

    app = build_app(walker)

    print(f"skywalker web UI: http://{host}:{port}/  (Ctrl-C to stop)")
    if host not in _LOOPBACK:
        print(f"  (remote access: ssh -L {port}:localhost:{port} "
             "user@this-host, then open http://localhost:%d/)" % port)

    if open_browser:
        import webbrowser
        threading.Timer(
            1.0, lambda: webbrowser.open(f"http://{host}:{port}/")).start()

    try:
        app.run(host=host, port=port, debug=debug)
    except OSError as exc:
        walker.logger.error(
            f"Could not start the web server on {host}:{port}: {exc}. "
            f"Try --web-port {port + 1}.")
        raise
