"""Interactive mouse-hover time cursor for the skywalker altitude plot.

Hovering the altitude panel drops a vertical line at the pointer time, slides a
circle along every plotted track at that time, and lists the time, altitude and
airmass of every object in a single annotation box. Hovering the polar skychart
snaps to the nearest plotted sample and drives the same readout, so the two
panels cross-reference each other.
"""

import numpy as np
import matplotlib

# Builtin matplotlib backends that cannot show an interactive window. Frozen as
# a literal on purpose: matplotlib.rcsetup.non_interactive_bk is deprecated
# since 3.9 and matplotlib.backends.backend_registry only exists from 3.9 on,
# while the package pins no matplotlib version.
NON_INTERACTIVE_BACKENDS = frozenset(
    ['agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template'])

# No value to report: outside the sampled range, or below the horizon.
NOVALUE = '--'


def hover_is_possible(fig):
    """Check whether the figure's backend can show an interactive window.

    Parameters:
    -----------
    fig : matplotlib.figure.Figure
        A figure that has already been created, so that the backend is
        resolved and querying it cannot trigger a GUI import.
    """
    if matplotlib.get_backend().lower() in NON_INTERACTIVE_BACKENDS:
        return False
    # Note: fig.canvas.manager is a real manager even under Agg, so it can
    # only ever be a secondary test.
    return getattr(fig.canvas, 'manager', None) is not None


def airmass(alt):
    """Airmass at altitude alt in degrees. Same formula as cli.set_plot()."""
    return 1. / np.cos(np.deg2rad(90. - alt))


class TimeCursor:
    """Vertical time cursor reporting altitude and airmass for all tracks.

    Parameters:
    -----------
    fig : matplotlib.figure.Figure
        The figure holding the plot.
    ax1 : matplotlib.axes.Axes
        The altitude panel. Hosts the vertical line, the sliding markers and
        the annotation box.
    tracks : list of dict
        One entry per plotted object, with keys 'name', 'alt', 'az', 'color'
        and 'is_moon'. 'alt' and 'az' are arrays in degrees, sampled on x.
    x : np.ndarray
        Hours from local midnight, strictly increasing, matching 'alt'/'az'.
    ax2 : matplotlib.axes.Axes, optional
        The airmass twin axis. Must be passed because it is created last and
        therefore receives nearly every mouse event over the panel.
    ax3 : matplotlib.axes.Axes, optional
        The polar skychart axes, or None when --skychart was not given.
    utcoffset : float, optional
        Offset of local time from UTC in hours, for the readout header.
    minalt : float, optional
        Minimum safe telescope altitude in degrees. Rows below it are flagged.
    moon_brightness : float, optional
        Moon illuminated fraction, 0-1. Appended as 'illum. NN%' to the
        Moon's own row when given. Default is None, which omits it.
    max_rows : int, optional
        Largest number of object rows to list. Default is 14.
    snap_px : float, optional
        How close the pointer must be, in pixels, to snap to a skychart
        sample. Default is 45.
    useblit : bool, optional
        Blit instead of redrawing the whole figure. Default is True.
    logger : logging.Logger, optional
        Used to report a blitting fallback.

    Notes:
    ------
    The clock time is derived from the pointer position rather than from the
    axis tick labels, which cli.set_plot() pins with a fixed formatter; after a
    toolbar zoom the ticks go stale while this readout stays correct.
    """

    def __init__(self, fig, ax1, tracks, x, ax2=None, ax3=None,
                 utcoffset=0., minalt=0., moon_brightness=None,
                 max_rows=14, snap_px=45., useblit=True, logger=None):

        self.fig = fig
        self.ax1 = ax1
        self.ax2 = ax2
        self.ax3 = ax3
        self.tracks = list(tracks)
        self.x = np.asarray(x, dtype=float)
        self.utcoffset = utcoffset
        self.minalt = minalt
        self.moon_brightness = moon_brightness
        self.max_rows = max_rows
        self.snap_px = snap_px
        self.useblit = bool(useblit) and fig.canvas.supports_blit
        self.logger = logger

        # Azimuth has to be unwrapped before it can be interpolated: a track
        # crossing 359 -> 1 deg would otherwise interpolate through 180 deg and
        # the skychart marker would jump across the chart.
        self._theta = [np.unwrap(np.deg2rad(_t['az'])) for _t in self.tracks]

        self._namew = min(12, max(len(_t['name']) for _t in self.tracks))
        self._cids = []
        self._background = None
        self._bg_wh = None
        self._pix = None
        self._visible = False

        self.set_artists()

    def set_artists(self):
        """Create the cursor artists, all invisible until the first hover."""
        self._vline = self.ax1.axvline(self.x[0], color='0.6', lw=0.8,
                                       zorder=20, visible=False,
                                       animated=self.useblit)
        self._marks = []
        self._polar_marks = []
        for _t in self.tracks:
            if _t['color'] is None:
                # The single-object track is coloured by the azimuth colormap,
                # so it has no one colour: a hollow ring reads over both
                # viridis and the near-black twilight fills.
                _style = {'mfc': 'none', 'mec': 'r', 'mew': 2.}
            else:
                _style = {'mfc': _t['color'], 'mec': 'k', 'mew': 0.8}
            self._marks.append(
                self.ax1.plot([], [], marker='o', ms=8, ls='none', zorder=21,
                              visible=False, animated=self.useblit,
                              **_style)[0])
            if self.ax3 is not None:
                self._polar_marks.append(
                    self.ax3.plot([], [], marker='o', ms=10, ls='none',
                                  mfc='none', mew=1.8, zorder=21,
                                  mec='r' if _t['color'] is None
                                  else _t['color'],
                                  visible=False, animated=self.useblit)[0])

        # A monospace font is load-bearing here, not cosmetic: the default
        # proportional face will not align the columns.
        self._ann = self.ax1.annotate(
            '', xy=(0., 0.), xytext=(14, 14), textcoords='offset points',
            ha='left', va='bottom', fontsize=8, family='monospace',
            color='k', zorder=22, visible=False, annotation_clip=False,
            animated=self.useblit,
            bbox=dict(boxstyle='round,pad=0.4', fc='w', ec='0.35',
                      alpha=0.92))

        self._artists = ([self._vline] + self._marks + self._polar_marks
                         + [self._ann])

    def connect(self):
        """Connect the cursor to the canvas and return self.

        The return value must be kept alive by the caller: matplotlib holds
        callbacks weakly, so an unreferenced cursor is garbage collected and
        the hover silently stops working.
        """
        _canvas = self.fig.canvas
        self._cids = [
            _canvas.mpl_connect('motion_notify_event', self.on_move),
            _canvas.mpl_connect('draw_event', self.set_background),
            _canvas.mpl_connect('axes_leave_event', self.on_leave),
            _canvas.mpl_connect('figure_leave_event', self.on_leave),
        ]
        return self

    def disconnect(self):
        """Disconnect every callback of this cursor."""
        for _cid in self._cids:
            self.fig.canvas.mpl_disconnect(_cid)
        self._cids = []

    def set_background(self, event=None):
        """Cache the clean figure for blitting. Bound to 'draw_event'."""
        self._pix = None            # the axes may have moved
        if not self.useblit or self.fig.canvas.is_saving():
            return
        self._background = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        self._bg_wh = self.fig.canvas.get_width_height()

    def on_leave(self, event):
        """Clear the cursor when the pointer leaves the axes or the figure."""
        self._hide()

    def on_move(self, event):
        """Update the cursor from a mouse motion event."""
        if not self.fig.canvas.widgetlock.available(self):
            return              # the toolbar is in pan or zoom mode
        _ax = event.inaxes
        if self.ax3 is not None and _ax is self.ax3:
            _xval, _yval = self._x_from_polar(event), None
        elif _ax is self.ax1 or (self.ax2 is not None and _ax is self.ax2):
            # ax2 is the airmass twin: it is created last, shares ax1's
            # position and wins event.inaxes over the whole panel. Reproject
            # the pixel through ax1 rather than trusting event.xdata, which is
            # only correct here because the two axes happen to share limits.
            _xval, _yval = self.ax1.transData.inverted().transform(
                (event.x, event.y))
        else:
            self._hide()        # the colorbar axes, or the figure margins
            return
        if _xval is None or not np.isfinite(_xval):
            self._hide()
            return
        self._update(_xval, yval=_yval)

    def _update(self, xval, yval=None):
        """Move every artist to xval and rebuild the readout."""
        _alts = [np.interp(xval, self.x, _t['alt'],
                           left=np.nan, right=np.nan) for _t in self.tracks]
        _keep = self._rows_to_keep(_alts)

        _rows = [self._clock_label(xval)]
        for _i, (_t, _alt) in enumerate(zip(self.tracks, _alts)):
            _has = np.isfinite(_alt) and _alt > 0.
            self._marks[_i].set_data([xval] if _has else [],
                                     [_alt] if _has else [])
            self._marks[_i].set_visible(_has)
            if self.ax3 is not None:
                _theta = np.interp(xval, self.x, self._theta[_i])
                self._polar_marks[_i].set_data(
                    [_theta % (2. * np.pi)] if _has else [],
                    [91. - _alt] if _has else [])
                self._polar_marks[_i].set_visible(_has)
            if _i in _keep:
                _rows.append(self._track_row(_t, _alt, _has))
        if len(_keep) < len(self.tracks):
            _rows.append('... +%i more' % (len(self.tracks) - len(_keep)))

        self._ann.set_text('\n'.join(_rows))
        self._vline.set_xdata([xval, xval])
        self._vline.set_visible(True)
        self._ann.set_visible(True)
        # An anchor is still needed when the pointer is on the skychart and
        # there is no y under it: mid-panel puts the box above and to the
        # right, deterministically.
        self._ann.xy = (xval, 45. if yval is None else yval)
        self._place_annotation()
        self._visible = True
        self._redraw()

    def _rows_to_keep(self, alts):
        """Pick which track rows fit in the box, keeping target-list order."""
        if len(self.tracks) <= self.max_rows:
            return set(range(len(self.tracks)))
        # Drop the rows with nothing to report first, then the lowest ones:
        # what is actually visible in the panel stays listed.
        _order = sorted(range(len(self.tracks)),
                        key=lambda i: (-alts[i] if np.isfinite(alts[i])
                                       and alts[i] > 0. else np.inf))
        return set(_order[:self.max_rows])

    def _track_row(self, track, alt, has_value):
        """Format one object row of the readout."""
        if not has_value:
            return '%-*s %5s %6s' % (self._namew, track['name'][:self._namew],
                                     NOVALUE, NOVALUE)
        _air = airmass(alt)
        # Blank out useless airmasses the same way the twin axis does.
        _airstr = '>10' if _air > 10. else '%.2f' % _air
        _flag = '!' if (not track['is_moon'] and alt < self.minalt) else ''
        _row = '%-*s %5.1f %6s%s' % (self._namew,
                                     track['name'][:self._namew],
                                     alt, _airstr, _flag)
        if track['is_moon'] and self.moon_brightness is not None:
            # Mirrors htmlplot's hovertemplate wording ('illum. NN%'), baked
            # in here as a scalar the same way: the fraction does not vary
            # with time, so it does not need a column of its own.
            _row += ' illum. %.0f%%' % (self.moon_brightness * 100.)
        return _row

    def _clock_label(self, xval):
        """Return the local clock time at xval, as the x ticks label it."""
        _hours = (xval + 24. if xval < 0. else xval) % 24.
        _hh = int(_hours)
        _mm = int(round((_hours - _hh) * 60.))
        if _mm == 60:                       # e.g. 01:59:40 reads as 02:00
            _mm = 0
            _hh = (_hh + 1) % 24
        return '%02d:%02d LT (UTC%+d)' % (_hh, _mm, int(self.utcoffset))

    def _place_annotation(self):
        """Flip and clamp the annotation box so it stays inside the panel."""
        _bb = self.ax1.bbox
        _xpix, _ypix = self.ax1.transData.transform(self._ann.xy)
        # Flip the offset and the alignment together, so the box always grows
        # away from the pointer and cannot hang off the edge it is nearest.
        _dx, _ha = (14, 'left') if (_xpix - _bb.x0) < 0.55 * _bb.width \
            else (-14, 'right')
        _dy, _va = (-12, 'top') if (_ypix - _bb.y0) > 0.5 * _bb.height \
            else (12, 'bottom')
        self._ann.set_ha(_ha)
        self._ann.set_va(_va)
        self._ann.xyann = (_dx, _dy)
        # The box grows with the number of rows, so it can still overflow
        # vertically after flipping. Measure and clamp; the extent is only
        # valid once the artist is visible and its text is set.
        _ext = self._ann.get_window_extent()
        _topoints = 72. / self.fig.dpi
        if _ext.y0 < _bb.y0:
            self._ann.xyann = (_dx, _dy + (_bb.y0 - _ext.y0) * _topoints)
        elif _ext.y1 > _bb.y1:
            self._ann.xyann = (_dx, _dy - (_ext.y1 - _bb.y1) * _topoints)

    def _polar_pixels(self):
        """Cache the display coords of every above-horizon skychart sample."""
        if self._pix is None:
            _pts, _idx = [], []
            for _i, _t in enumerate(self.tracks):
                _mask = _t['alt'] > 0.
                if not _mask.any():
                    continue
                _pts.append(self.ax3.transData.transform(np.column_stack(
                    [self._theta[_i][_mask] % (2. * np.pi),
                     91. - _t['alt'][_mask]])))
                _idx.append(np.flatnonzero(_mask))
            self._pix = ((np.vstack(_pts), np.concatenate(_idx))
                         if _pts else (None, None))
        return self._pix

    def _x_from_polar(self, event):
        """Return the time of the skychart sample nearest to the pointer.

        Nearest is measured in display pixels, not in polar data coordinates:
        the radial axis is inverted and angular distance scales with radius,
        so a fixed data distance would mean very different visual distances
        near the centre and near the rim.
        """
        _pts, _idx = self._polar_pixels()
        if _pts is None:
            return None
        _dist = ((_pts[:, 0] - event.x) ** 2 + (_pts[:, 1] - event.y) ** 2)
        _near = int(_dist.argmin())
        if _dist[_near] > self.snap_px ** 2:
            return None                 # not near any track
        return self.x[_idx[_near]]

    def _hide(self):
        """Hide every artist and commit the clean figure."""
        if not self._visible:
            return
        for _artist in self._artists:
            _artist.set_visible(False)
        self._visible = False
        # A blit would leave the last cursor painted on the real canvas.
        self.fig.canvas.draw_idle()

    def _redraw(self):
        """Blit the cursor over the cached background, or redraw."""
        _canvas = self.fig.canvas
        if (self.useblit and self._background is not None
                and _canvas.get_width_height() == self._bg_wh):
            try:
                _canvas.restore_region(self._background)
                for _artist in self._artists:
                    _artist.axes.draw_artist(_artist)
                _canvas.blit(self.fig.bbox)
                return
            except Exception as e:
                # Raising inside a callback would print a traceback on every
                # mouse motion, which is far worse than a slow cursor.
                if self.logger is not None:
                    self.logger.warning(
                        f"Hover blitting failed ({e}): falling back to full \
                        redraws.")
                self.useblit = False
                for _artist in self._artists:
                    _artist.set_animated(False)
        _canvas.draw_idle()
