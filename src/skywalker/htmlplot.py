"""Interactive HTML renderer for the skywalker altitude plot, using plotly.

Consumes the same track registry that cli.Skywalker.set_plot() builds while
plotting the matplotlib figure, so no astronomy is recomputed here: this
module only reshapes already-computed arrays into plotly traces.

plotly is an optional dependency: importing this module never requires it,
only calling render()/build_figure() does, and the ImportError raised then
points the user at how to install it.
"""

from datetime import timedelta

import numpy as np

from .plotdata import airmass_from_alt, mask_to_intervals, PALETTE as _PALETTE

_COMPASS_DEG = list(range(0, 360, 45))
_COMPASS_LABELS = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']

# matplotlib single-letter colour codes are not valid CSS colours, unlike the
# hex/named colours everywhere else in this module.
_MPL_COLOR = {'c': 'cyan', 'k': 'black', 'm': 'magenta',
             'r': 'red', 'g': 'green', 'b': 'blue', 'w': 'white'}


def _css_color(color):
    """Map a matplotlib single-letter colour code to a CSS colour name."""
    return _MPL_COLOR.get(color, color)


def _require_plotly():
    """Import plotly lazily, raising a clear error if it is not installed."""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as exc:
        raise ImportError(
            "--savehtml requires plotly, which is not installed.\n"
            "  Install it with:  pip install 'plotly>=5.15'\n"
            "  Or reinstall skywalker with the html extra:  "
            "pip install -e '.[html]'") from exc
    return go, make_subplots


def _track_color(track, index):
    """The track's own colour, or the next one from the fallback palette."""
    return _css_color(track['color']) if track['color'] \
        else _PALETTE[index % len(_PALETTE)]


def _local_midnight(local_times, delta_hours):
    """Recover local midnight from one sample of the shared time grid."""
    return local_times[0] - timedelta(hours=float(delta_hours[0]))


def build_figure(tracks, local_times, delta_hours, sun_alt, moon_alt,
                 moon_brightness, minalt, utcoffset_h, sitename,
                 nightstarts, make_skychart=False):
    """Build the interactive plotly figure.

    Parameters:
    -----------
    tracks : list of dict
        The same per-object registry cli.Skywalker.set_plot() builds: one
        entry per plotted object plus the Moon, each with 'name', 'alt',
        'az', 'color', 'is_moon', 'has_block', 'block_starts', 'block_ends',
        'moon_distance', 'label_x', 'label_y', 'label_color', 'chart_alt',
        'chart_az' and 'chart_hours' (the last three only when
        make_skychart is True).
    local_times : sequence of datetime.datetime
        Local wall-clock time of every sample, aligned with 'alt'/'az' and
        with sun_alt/moon_alt (500 samples).
    delta_hours : np.ndarray
        Decimal hours from local midnight, aligned with local_times. Used
        only to slice the observing-block spans; the axis itself uses
        local_times.
    sun_alt, moon_alt : np.ndarray
        Sun and Moon altitude in degrees, aligned with local_times.
    moon_brightness : float
        Moon illumination fraction, 0-1.
    minalt : float
        Minimum safe telescope altitude in degrees.
    utcoffset_h : int
        Local time's offset from UTC, in hours.
    sitename : str
    nightstarts : str
    make_skychart : bool, optional
        Whether to add the polar skychart panel. Default is False.
    """
    go, make_subplots = _require_plotly()

    local_times = np.asarray(local_times)
    delta_hours = np.asarray(delta_hours, dtype=float)
    sun_alt = np.asarray(sun_alt, dtype=float)
    moon_alt = np.asarray(moon_alt, dtype=float)
    _midnight = _local_midnight(local_times, delta_hours)

    if make_skychart:
        fig = make_subplots(
            rows=1, cols=2,
            specs=[[{'type': 'xy'}, {'type': 'polar'}]],
            column_widths=[0.62, 0.38])
    else:
        fig = go.Figure()

    # Row/col kwargs for the xy panel: only meaningful once make_subplots()
    # has actually been used, i.e. when the skychart panel exists too.
    _rc = {'row': 1, 'col': 1} if make_skychart else {}

    def _add(trace, col=1):
        if make_skychart:
            fig.add_trace(trace, row=1, col=col)
        else:
            fig.add_trace(trace)

    # Twilight/darkness bands, mirroring cli.Skywalker.set_plot() exactly.
    _bands = [
        ((sun_alt < 0.) & (sun_alt > -6.3), 'indigo', 0.8),
        ((sun_alt < -6.) & (sun_alt > -12.3), 'indigo', 0.9),
        ((sun_alt < -12.) & (sun_alt > -18.), 'indigo', 1.0),
        ((moon_alt < 0.) & (sun_alt < -18.), 'black', 1.0),
        ((moon_alt > 0.) & (sun_alt < -18.), 'midnightblue',
         1. - moon_brightness),
    ]
    for _mask, _color, _alpha in _bands:
        for _t0, _t1 in mask_to_intervals(local_times, _mask):
            fig.add_shape(type='rect', xref='x', yref='y domain',
                         x0=_t0, x1=_t1, y0=0, y1=1,
                         fillcolor=_color, opacity=_alpha, line_width=0,
                         layer='below', **_rc)

    # minalt limit line.
    if minalt > 1:
        fig.add_hline(y=minalt, line=dict(color='red', dash='dash'), **_rc)

    _obj_index = 0
    for _track in tracks:
        _alt = np.asarray(_track['alt'], dtype=float)
        _az = np.asarray(_track['az'], dtype=float)
        _air = airmass_from_alt(_alt)
        _y = np.where(_alt > 0., _alt, np.nan)
        _customdata = np.stack([_air, _az], axis=-1)

        _name = _track['name']
        if _track['is_moon']:
            _color = _css_color(_track['color'])
            # No <extra></extra>: this trace is on the same unified x-hover
            # panel as the objects, so trace.name must still supply the row
            # label (<extra> would replace it, see htmlplot module notes).
            # Built with '+', not '%', since the template itself already
            # uses a '%{...}' syntax that a Python '%' operator would choke
            # on (the moon illumination is baked in as a literal, since it
            # is a single scalar rather than a per-sample value).
            _hover = ('%{y:.1f}&deg;  &middot;  illum. '
                     + ('%.0f' % (moon_brightness * 100.)) + '%')
        else:
            _color = _track_color(_track, _obj_index)
            _obj_index += 1
            _hover = ('%{y:.1f}&deg;  &middot;  X %{customdata[0]:.2f}'
                      '  &middot;  Az %{customdata[1]:.0f}&deg;')

        _add(go.Scatter(x=local_times, y=_y, mode='lines', name=_name,
                        legendgroup=_name,
                        line=dict(color=_color,
                                  dash='dash' if _track['is_moon'] else
                                  'solid'),
                        customdata=_customdata, hovertemplate=_hover))

        if _track['has_block']:
            _mask = ((delta_hours >= _track['block_starts'])
                     & (delta_hours <= _track['block_ends']))
            if _mask.any():
                _add(go.Scatter(x=local_times[_mask], y=_y[_mask],
                               mode='lines', fill='tozeroy',
                               fillcolor=_color, opacity=0.35,
                               line_width=0, showlegend=False,
                               name=_name, legendgroup=_name,
                               hoverinfo='skip'))

        if _track['moon_distance'] is not None:
            fig.add_annotation(
                x=_midnight + timedelta(hours=float(_track['label_x'])),
                y=_track['label_y'],
                text='%i' % _track['moon_distance'], showarrow=False,
                font=dict(size=10, color=_css_color(_track['label_color'])),
                name=_track['name'], **_rc)

    if make_skychart:
        _obj_index = 0
        for _track in tracks:
            # Mirror the colour-assignment loop above exactly, so a track's
            # skychart marker always matches its altitude-panel line.
            if not _track['is_moon']:
                _color = _track_color(_track, _obj_index)
                _obj_index += 1
            else:
                _color = _css_color(_track['color'])
            if _track.get('chart_alt') is None:
                continue
            _calt = np.asarray(_track['chart_alt'], dtype=float)
            _caz = np.asarray(_track['chart_az'], dtype=float)
            _chours = np.asarray(_track['chart_hours'], dtype=float)
            fig.add_trace(go.Scatterpolar(
                r=91. - _calt, theta=_caz, mode='markers+text',
                text=['%.1f' % h for h in _chours],
                textposition='top center', textfont=dict(size=8),
                marker=dict(size=9, symbol='star', color=_color),
                name=_track['name'], legendgroup=_track['name'],
                showlegend=False,
                customdata=_calt,
                hovertemplate=('<b>%{fullData.name}</b><br>'
                               'Alt %{customdata:.1f}&deg;  '
                               'Az %{theta:.0f}&deg;<extra></extra>')),
                row=1, col=2)

    _title = f"Night starts: {nightstarts} @ {sitename}"
    fig.update_layout(
        title=_title,
        uirevision='skywalker',       # keep zoom/pan across figure rebuilds
        hovermode='x unified',
        hoverlabel=dict(font_family='monospace', align='left',
                        namelength=-1),
        legend=dict(x=0.99, y=0.99, xanchor='right', yanchor='top'),
        margin=dict(l=60, r=40, t=60, b=50),
        width=1400 if make_skychart else 900, height=620)

    fig.update_xaxes(
        title=f'Local Time [UTC{utcoffset_h:+d}]',
        showspikes=True, spikemode='across', spikesnap='cursor',
        spikethickness=1, spikedash='dot', spikecolor='#555',
        hoverformat='%H:%M', **_rc)
    fig.update_yaxes(title='Altitude [deg]', range=[0, 90], **_rc)

    if make_skychart:
        fig.update_polars(
            angularaxis=dict(direction='counterclockwise', rotation=90,
                             tickmode='array', tickvals=_COMPASS_DEG,
                             ticktext=_COMPASS_LABELS),
            radialaxis=dict(range=[1, 91], angle=-45))

    return fig


def write_html(fig, filename, include_js='embed'):
    """Write fig to filename. Returns the path written.

    Parameters:
    -----------
    fig : plotly.graph_objects.Figure
    filename : str
    include_js : str, optional
        'embed' (self-contained, works offline), 'cdn' (small file, needs
        internet) or 'directory' (a shared plotly.min.js file alongside the
        output). Default is 'embed'.
    """
    _include = {'embed': True, 'cdn': 'cdn', 'directory': 'directory'}[include_js]
    fig.write_html(filename, include_plotlyjs=_include, auto_open=False)
    return filename


def render(tracks, local_times, delta_hours, sun_alt, moon_alt,
          moon_brightness, minalt, utcoffset_h, sitename, nightstarts,
          filename, make_skychart=False, include_js='embed', logger=None):
    """Build the figure and write it to filename. Returns the path written."""
    fig = build_figure(tracks, local_times, delta_hours, sun_alt, moon_alt,
                       moon_brightness, minalt, utcoffset_h, sitename,
                       nightstarts, make_skychart=make_skychart)
    return write_html(fig, filename, include_js=include_js)
