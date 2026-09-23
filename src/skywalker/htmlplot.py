"""Interactive HTML renderer for the skywalker altitude plot, using plotly.

Consumes the same track registry that cli.Skywalker.set_plot() builds while
plotting the matplotlib figure, so no astronomy is recomputed here: this
module only reshapes already-computed arrays into plotly traces.

plotly is an optional dependency: importing this module never requires it,
only calling render()/build_figure() does, and the ImportError raised then
points the user at how to install it.
"""

from datetime import datetime, timedelta

import numpy as np

from .plotdata import (airmass_from_alt, astro_night_mask, mask_to_intervals,
                       PALETTE as _PALETTE)

_COMPASS_DEG = list(range(0, 360, 45))
_COMPASS_LABELS = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']

# matplotlib single-letter colour codes are not valid CSS colours, unlike the
# hex/named colours everywhere else in this module.
_MPL_COLOR = {'c': 'cyan', 'k': 'black', 'm': 'magenta',
             'r': 'red', 'g': 'green', 'b': 'blue', 'w': 'white'}


def _css_color(color):
    """Map a matplotlib single-letter colour code to a CSS colour name."""
    return _MPL_COLOR.get(color, color)


def _over_white(color, alpha):
    """Opaque (r, g, b) of color laid at alpha over white, 0-255 ints.

    The night bands and the skychart are pre-blended rather than drawn
    translucent, so the page theme's plot background never shows through
    them: a bright-Moon night reads pale and a dark one deep blue in both
    the light and the dark theme, as on the matplotlib figure's white axes.
    """
    from matplotlib.colors import to_rgb
    return tuple(int(round(255. * (alpha * _c + 1. - alpha)))
                 for _c in to_rgb(color))


def _rgb_css(rgb):
    """CSS 'rgb(r, g, b)' string for an (r, g, b) tuple of 0-255 ints."""
    return 'rgb(%d, %d, %d)' % rgb


def _is_dark(rgb):
    """Whether light text/gridlines are needed on this (r, g, b) fill."""
    _r, _g, _b = rgb
    return (0.2126 * _r + 0.7152 * _g + 0.0722 * _b) / 255. < 0.5


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
                 nightstarts, make_skychart=False, dark=False):
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
    dark : bool, optional
        Render for a dark browser theme: switches the plotly template to
        'plotly_dark', fixes paper/plot background to '#1e1e1e', and
        lightens the hover label and axis-spike colours. The night bands
        and the skychart background do not follow the theme: they are
        coloured by the Sun and the Moon's illumination alone. Default is
        False, which is what render()/write_html() and every --savehtml
        path use.
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

    # Twilight/darkness bands, mirroring cli.Skywalker.set_plot()'s colours
    # and alphas, but pre-blended over white (see _over_white()) and drawn
    # opaque, so they look the same whatever the page theme. The Moon-up
    # astronomical night is therefore as pale as the Moon is bright.
    _bands = [
        ((sun_alt < 0.) & (sun_alt > -6.3), 'indigo', 0.8),
        ((sun_alt < -6.) & (sun_alt > -12.3), 'indigo', 0.9),
        ((sun_alt < -12.) & (sun_alt > -18.), 'indigo', 1.0),
        ((moon_alt < 0.) & astro_night_mask(sun_alt), 'black', 1.0),
        ((moon_alt > 0.) & astro_night_mask(sun_alt), 'midnightblue',
         1. - moon_brightness),
    ]
    for _mask, _color, _alpha in _bands:
        for _t0, _t1 in mask_to_intervals(local_times, _mask):
            fig.add_shape(type='rect', xref='x', yref='y domain',
                         x0=_t0, x1=_t1, y0=0, y1=1,
                         fillcolor=_rgb_css(_over_white(_color, _alpha)),
                         opacity=1., line_width=0,
                         layer='below', **_rc)

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
            _hover = ('%{y:.1f}&deg;  ·  illum. '
                     + ('%.0f' % (moon_brightness * 100.)) + '%')
        else:
            _color = _track_color(_track, _obj_index)
            _obj_index += 1
            _hover = ('%{y:.1f}&deg;  ·  X %{customdata[0]:.2f}'
                      '  ·  Az %{customdata[1]:.0f}&deg;')

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

    # minalt limit line. Added only now, after the track loop above has
    # added at least one go.Scatter trace to this subplot: plotly 7's
    # add_hline(row=, col=) silently fails to append any shape at all --
    # no exception, just a no-op -- when called against a make_subplots()
    # panel that does not yet hold a trace (add_shape() has no such
    # requirement, which is why the twilight/darkness bands above, added
    # before any trace exists, render fine). This was a real bug: with
    # make_skychart=True the red minalt line never appeared, silently,
    # because it used to be added here before the loop. tracks is never
    # empty in practice (ordered_tracks() always includes the Moon unless
    # the caller explicitly filters it out), but a caller that manages to
    # pass an empty tracks list will still lose the line -- not fixed
    # here, since nothing in this codebase currently does that.
    if minalt > 1:
        fig.add_hline(y=minalt, line=dict(color='red', dash='dash'), **_rc)

    if make_skychart:
        # The whole chart is tinted by the Moon's illuminated fraction:
        # black at alpha 1 - brightness over white, as the matplotlib
        # skychart draws it, so a new Moon gives a black sky and a full
        # one a white sky. Like the staralt bands it does not follow the
        # page theme; labels and gridlines pick light or dark ink to stay
        # readable on it.
        _sky = _over_white('black', 1. - moon_brightness)
        _ink = '#e0e0e0' if _is_dark(_sky) else '#222222'
        _grid = ('rgba(255, 255, 255, 0.3)' if _is_dark(_sky)
                 else 'rgba(0, 0, 0, 0.2)')
        # Red ring below minalt, as the matplotlib skychart draws it
        # (red at alpha 0.7), pre-blended over white like the sky itself.
        # Added before the stars so they draw on top of it.
        if minalt > 0.:
            fig.add_trace(go.Barpolar(
                r=[minalt], base=[91. - minalt], theta=[0.], width=[360.],
                marker=dict(color=_rgb_css(_over_white('red', 0.7)),
                            line_width=0),
                showlegend=False, hoverinfo='skip', name='minalt',
                meta='sw-fixed'),
                row=1, col=2)
        _obj_index = 0
        for _track in tracks:
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
                textposition='top center',
                textfont=dict(size=8, color=_ink),
                marker=dict(size=9, symbol='star', color=_color),
                name=_track['name'], legendgroup=_track['name'],
                showlegend=False,
                customdata=_calt,
                hovertemplate=('<b>%{fullData.name}</b><br>'
                               'Alt %{customdata:.1f}&deg;  '
                               'Az %{theta:.0f}&deg;<extra></extra>')),
                row=1, col=2)

    _title = f"Night starts: {nightstarts} @ {sitename}"
    _hoverlabel = dict(font_family='monospace', align='left',
                       namelength=-1)
    _layout_kwargs = dict(
        title=_title,
        uirevision='skywalker',       # keep zoom/pan across figure rebuilds
        hovermode='x unified',
        hoverlabel=_hoverlabel,
        legend=dict(x=0.99, y=0.99, xanchor='right', yanchor='top'),
        margin=dict(l=60, r=40, t=60, b=50),
        width=1400 if make_skychart else 900, height=620)
    if dark:
        # Not folded into _hoverlabel/_layout_kwargs above: keeping the
        # light-mode call's kwargs untouched means dark=False (every
        # --savehtml/render()/write_html() caller) reaches
        # fig.update_layout() with the exact same arguments as before
        # this parameter existed.
        _hoverlabel.update(bgcolor='#2a2a2a', font=dict(color='#e0e0e0'))
        _layout_kwargs.update(template='plotly_dark',
                              paper_bgcolor='#1e1e1e',
                              plot_bgcolor='#1e1e1e')
    fig.update_layout(**_layout_kwargs)

    fig.update_xaxes(
        title=f'Local Time [UTC{utcoffset_h:+d}]',
        showspikes=True, spikemode='across', spikesnap='cursor',
        spikethickness=1, spikedash='dot',
        spikecolor='#bbb' if dark else '#555',
        hoverformat='%H:%M', **_rc)
    fig.update_yaxes(title='Altitude [deg]', range=[0, 90], **_rc)

    if make_skychart:
        fig.update_polars(
            bgcolor=_rgb_css(_sky),
            angularaxis=dict(direction='counterclockwise', rotation=90,
                             tickmode='array', tickvals=_COMPASS_DEG,
                             ticktext=_COMPASS_LABELS, gridcolor=_grid),
            radialaxis=dict(range=[1, 91], angle=-45, gridcolor=_grid,
                            tickfont=dict(color=_ink)))

    return fig


def build_year_figure(curves, dates, minalt, sitename, marked_date,
                      metric='alt', dark=False):
    """Build the year-view plotly figure: one line per target, per month.

    Pure reshaper, like build_figure(): it does no astronomy, only lays out
    the per-month peak altitudes and usable hours that
    cli.Skywalker.year_max_altitudes() has already computed, one line per
    target.

    Parameters:
    -----------
    curves : list of dict
        One entry per target, each with 'name', 'color' (str or None, same
        convention as a track's 'color'), 'peak_alt' (np.ndarray, shape
        (n,), NaN for a month where the target is not observable during
        astronomical night), 'peak_time' (list of str, length n, 'HH:MM' or
        '—' where peak_alt is NaN) and 'hours_up' (np.ndarray, shape (n,)).
    dates : list of datetime.date
        One reference date per month, length n, aligned with every curve's
        arrays.
    minalt : float
        Minimum safe telescope altitude in degrees.
    sitename : str
    marked_date : datetime.date
        The year-view date box's chosen date. Drawn as a vertical dotted
        line across the figure and labelled with its ISO date, replacing
        the old fixed line at dates[0] ("tonight"). May fall on any day,
        not just a date already present in dates.
    metric : str, optional
        'alt' (default, and the only behaviour before this parameter
        existed) plots peak_alt on the y axis: minalt is drawn as a red
        dashed hline and the y axis is fixed to [0, 90]. 'hours' plots
        hours_up instead: the minalt hline is dropped -- minalt is the
        threshold hours_up was counted against, not a y value -- and the y
        axis autoscales from zero (rangemode='tozero') instead of being
        fixed, since the natural ceiling is the length of the night and
        varies by site and season. Either way the hover shows both
        peak_alt and hours_up; only which one drives the y axis, and which
        one sits in customdata[1], swaps. The y axis's own uirevision (see
        below) is keyed by metric as well as year, so switching between
        'alt' and 'hours' always resets the y axis to its default
        range/autorange instead of keeping a zoom taken under the other
        metric's incompatible scale.
    dark : bool, optional
        Render for a dark browser theme: 'plotly_dark' template, fixed
        '#1e1e1e' paper/plot background, a lightened hover label, and a
        lighter '#bbb' marked-date line (the light-mode '#888' reads as
        near-invisible against a dark paper). Default is False, which
        reproduces today's figure exactly.

    The y value in 'alt' mode is the peak altitude reached during
    astronomical night (Sun below -18 deg, see plotdata.astro_night_mask()),
    which is a strictly tighter cut than the 0 deg "Sun below the horizon"
    one behind the "Peak alt" column of the web UI's table (track_summary()'s
    night_mask). The first (tonight's) point on this curve can therefore
    legitimately read slightly lower than that column's value.

    Gap semantics differ between the two metrics -- this is deliberate, not
    a bug to "fix" later. In 'alt' mode a non-observable month is NaN in
    peak_alt and renders as a real gap (connectgaps=False keeps the line
    from bridging it). In 'hours' mode that same month's hours_up is 0.0,
    not NaN: a target that is up for zero usable hours that month is a
    meaningful, honest zero, so it is plotted as a real point at y=0 and
    must never be converted to a gap. The "never observable all year"
    legend treatment (name suffix and visible='legendonly') keys off
    peak_alt being all-NaN in both modes, so a target's legend entry does
    not change meaning when the metric is flipped.
    """
    go, _ = _require_plotly()

    fig = go.Figure()

    _obj_index = 0
    for _curve in curves:
        _peak_alt = np.asarray(_curve['peak_alt'], dtype=float)
        _hours_up = np.asarray(_curve['hours_up'], dtype=float)
        _peak_time = _curve['peak_time']
        _color = _track_color(_curve, _obj_index)
        _obj_index += 1

        _name = _curve['name']
        _visible = True
        if np.all(np.isnan(_peak_alt)):
            _name = _name + ' (never observable)'
            _visible = 'legendonly'

        if metric == 'hours':
            _y = _hours_up
            _customdata = np.stack([_peak_time, _peak_alt], axis=-1)
            _hover = ('%{y:.1f} h usable · peak '
                     '%{customdata[1]:.1f}&deg; at %{customdata[0]}')
        else:
            _y = _peak_alt
            _customdata = np.stack([_peak_time, _hours_up], axis=-1)
            _hover = ('%{y:.1f}&deg; at %{customdata[0]} · '
                     '%{customdata[1]:.1f} h usable')

        fig.add_trace(go.Scatter(
            x=dates, y=_y, mode='lines+markers', name=_name,
            legendgroup=_name, connectgaps=False, visible=_visible,
            line=dict(color=_color), marker=dict(color=_color),
            customdata=_customdata, hovertemplate=_hover))

    if metric != 'hours' and minalt > 1:
        fig.add_hline(y=minalt, line=dict(color='red', dash='dash'))

    # add_vline raises on a bare datetime.date on this plotly version (it
    # tries to add an int offset to it internally); a full datetime is
    # accepted, so promote marked_date before handing it over, and fall
    # back to an equivalent add_shape()/add_annotation() pair if some
    # other installed version still rejects it.
    _marked_dt = datetime.combine(marked_date, datetime.min.time())
    _marked_label = marked_date.strftime('%Y-%m-%d')
    _marker_color = '#bbb' if dark else '#888'
    try:
        fig.add_vline(x=_marked_dt,
                     line=dict(color=_marker_color, dash='dot'),
                     annotation_text=_marked_label)
    except (TypeError, ValueError):
        fig.add_shape(type='line', xref='x', yref='y domain',
                     x0=_marked_dt, x1=_marked_dt, y0=0, y1=1,
                     line=dict(color=_marker_color, dash='dot'))
        fig.add_annotation(x=_marked_dt, y=1, yref='y domain',
                          yanchor='bottom', text=_marked_label, showarrow=False)

    if metric == 'hours':
        _title = ("Year view: hours usable during astronomical night "
                 f"@ {sitename}")
    else:
        _title = ("Year view: peak altitude during astronomical night "
                 f"@ {sitename}")
    _hoverlabel = dict(font_family='monospace', align='left',
                       namelength=-1)
    _layout_kwargs = dict(
        title=_title,
        # Keyed by year, not a bare constant: this governs the x axis,
        # the legend and per-trace visibility (legendonly toggles) --
        # a target-selection or marked-date change keeps the same key,
        # so those survive such a refresh, but a year switch changes
        # the key and resets them, since the old state is almost
        # certainly meaningless against a different year's data. The y
        # axis does NOT use this key -- see its own uirevision below,
        # which also resets on a metric switch; 'alt' and 'hours' have
        # incompatible y scales, so replaying one metric's zoom onto
        # the other silently produced a nonsensical axis. That was the
        # bug: before this uirevision was added, the y axis inherited
        # this same year-only key and kept a stale zoomed range across
        # an alt<->hours switch.
        uirevision=f'skywalker-year-{dates[0].year}',
        hovermode='x unified',
        hoverlabel=_hoverlabel,
        height=360, autosize=True, width=None,
        margin=dict(l=60, r=40, t=50, b=40))
    if dark:
        _hoverlabel.update(bgcolor='#2a2a2a', font=dict(color='#e0e0e0'))
        _layout_kwargs.update(template='plotly_dark',
                              paper_bgcolor='#1e1e1e',
                              plot_bgcolor='#1e1e1e')
    fig.update_layout(**_layout_kwargs)

    # The y axis gets its own uirevision, keyed by year AND metric,
    # separate from the figure-wide one above (year only): 'alt' and
    # 'hours' have incompatible y scales (0-90 deg fixed range vs. an
    # autoranged, site/season-dependent hour count), so a zoom taken in
    # one metric must not be replayed onto the other's axis when the
    # metric switch rebuilds this figure with the same top-level
    # uirevision. Keeping it on the y axis only (not the shared
    # uirevision, not the x axis) means a target-selection or
    # marked-date change still keeps the user's zoom within one metric,
    # exactly as before -- only crossing the alt<->hours boundary now
    # forces the axis back to its default range/autorange.
    _y_uirevision = f'skywalker-year-{dates[0].year}-{metric}'
    fig.update_xaxes(title='Date', tickformat='%b %Y', dtick='M1')
    if metric == 'hours':
        fig.update_yaxes(
            title=f'Hours above {minalt:.0f} deg in astro. night',
            rangemode='tozero', uirevision=_y_uirevision)
    else:
        fig.update_yaxes(title='Peak altitude in astro. night [deg]',
                         range=[0, 90], uirevision=_y_uirevision)

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
