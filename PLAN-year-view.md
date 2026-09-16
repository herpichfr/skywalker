# Year planner: seasonal visibility panel in the web UI

## Context

`--web` currently answers one question: *what does tonight look like?* It shows a
staralt panel plus an optional skychart (`htmlplot.build_figure()`), and a target
table underneath. There is no way to ask the other question an observer always has:
*which month should I propose this target for?* Today that means re-running
`skywalker` twelve times with twelve different `-ns` dates and eyeballing the PNGs.

This adds a third panel, directly under the existing plot, that answers it in one
look: date on the x-axis spanning a full year from the current `--nightstarts`,
altitude on the y-axis, one line per target, and each point is that target's
**highest altitude reached while the Sun is below -18 deg** on that night. A month
where the target never clears the horizon during astronomical night is a gap in the
line, not a zero — so the shape of the curve *is* the observing season.

Decisions already taken: 13 monthly samples (start date plus twelve months, so the
year is fully spanned); the panel sits directly under the main plot and above the
add-target controls; and it is **lazy** — nothing is computed until the user clicks
a toggle, so app startup stays as fast as it is now.

## Design in one paragraph

All the astronomy is one broadcast `transform_to()` per target. Build a single flat
`AltAz` frame of 13 nights x 145 samples (10-minute cadence over a 24 h window
centred on each night's local midnight), evaluate the Sun on it **once** to get the
astronomical-night mask, then for each target do one `SkyCoord.transform_to()` over
the whole flat grid, reshape to `(13, 145)`, mask to night, and `nanmax` along the
sample axis. Cache the 13-element result per target in the registry. Rendering is a
new plotly figure builder next to the existing one; the Dash side is two new
callbacks that obey the module's existing acyclic mutator/consumer split.

## Files to modify

### 0. `src/skywalker/plotdata.py` — one shared definition of astronomical night

Right now `-18` appears exactly once in the repo, inlined in the twilight-band list
of `htmlplot.build_figure()` (htmlplot.py:122-125). There is no reusable helper, and
nothing else uses the astronomical-night threshold — `TrackRegistry.night_mask`
(webapp.py:147) and `compute_track()`'s observability mask (cli.py:557-558) both cut
at the Sun below the *horizon*, 0 deg, which is a different and looser thing.

Add to `plotdata.py` (which is pure numpy by design, so this adds no import):

```python
ASTRO_NIGHT_ALT_DEG = -18.0

def astro_night_mask(sun_alt):
    """True where the Sun is below the astronomical-twilight limit."""
```

and refactor the two `sun_alt < -18.` literals in `htmlplot.py`'s `_bands` to call
it, so the year panel and the dark bands of the main panel can never drift apart.

### 1. `src/skywalker/cli.py` — the astronomy (new methods on `Skywalker`)

This is the only module that imports astropy/astroplan for computation, so the new
math belongs here, mirroring `set_night_frames()` (cli.py:344) and `compute_track()`
(cli.py:525).

**`set_year_frames(self, n_months=13, samples=145, sun_alt_limit=-18.)`**

- Step the night start dates with `pd.Timestamp(self.nightstarts) +
  pd.DateOffset(months=i)` for `i in range(n_months)`. `pandas` is already imported
  (cli.py:13) and `DateOffset` clamps month-end correctly (Jan 31 -> Feb 28).
- For each date reuse the exact local-midnight convention of `set_night_frames()`:
  `_night_ends = (Time(date + "T" + self.inithour) - self.utcoffset + .5*u.day)`,
  then `_midnight = Time(f"{_night_ends:%Y-%m-%d}T00:00:00") - self.utcoffset`. This
  guarantees sample 0 of the year grid is the same night the main panel is drawing.
- Concatenate all `n_months * samples` epochs into **one** `Time` array, build
  **one** `AltAz(obstime=..., location=self.location)`, and call `get_body('sun',
  ...)` **once** on it. Reshape the resulting altitudes to `(n_months, samples)`.
- Store: `self.year_dates` (list of `datetime.date`), `self.year_frame` (the flat
  `AltAz`), `self.year_shape`, `self.year_night_mask` (bool `(n_months, samples)`,
  from `plotdata.astro_night_mask()`), `self.year_local_times` (`(n_months, samples)` local
  wall-clock, for hover labels).
- Idempotent: return immediately if `self.year_frame` is already set.
- Acceptance: `year_night_mask.sum(axis=1)` is > 0 for every month at a mid-latitude
  site such as the OPD example, and month 0's mask matches the astronomical-night
  band of the existing `sunaltaz_time_overnight` to within one grid step.

**`year_max_altitudes(self, obj_coords)`**

- Takes a `SkyCoord` that may hold **one or many** targets, so the caller can batch —
  one code path for "app just opened with 30 targets" and "user added one more".
- `set_year_frames()` first (self-guarding), then a single broadcast transform:
  reshape the coordinates to `(N, 1)` and `transform_to(self.year_frame)` against the
  flat `(n_months*samples,)` obstime, giving `(N, n_months*samples)` from **one**
  ERFA-vectorised call; reshape to `(N, n_months, samples)`.
- `np.where(self.year_night_mask, alt, np.nan)` (broadcasting the `(n_months,
  samples)` mask over the target axis), then `np.nanmax(..., axis=2)` under
  `warnings.catch_warnings()` (an all-NaN month is expected and must not warn), then
  `np.where(peak > 0., peak, np.nan)` so a target below the horizon all night is a
  gap, not a point at 0.
- Also return, per target per month: `peak_time` (local `HH:MM` at the argmax, `'—'`
  when NaN) and `hours_up` (`night_mask & (alt > self.minalt)` count x grid step) —
  both feed the hover box, and `hours_up` is what actually tells the user how much of
  the month is usable, at no extra cost.
- Returns a list of `N` dicts keyed `peak_alt`, `peak_time`, `hours_up`.
- Acceptance: for the example file's NGC253 at OPD, `peak_alt` peaks in the
  September-November range and goes NaN or low around April-May; month 0's value
  equals the `peakalt` the table already shows for that night (same definition,
  modulo the -18 deg cut vs. the table's sun-below-horizon cut — see Risks).

### 2. `src/skywalker/htmlplot.py` — the renderer

**`build_year_figure(curves, dates, minalt, sitename, nightstarts)`**

Goes in `htmlplot.py` rather than a new `yearplot.py`: it reuses `_require_plotly()`
(htmlplot.py:32), `_css_color()` and `_track_color()` (htmlplot.py:46) verbatim, and
a separate module would have to import them or duplicate them. `hover.py` is *not*
involved — it is matplotlib-only cursor machinery used solely by the non-web GUI
path (`Skywalker.set_hover()`, cli.py:868); the web UI's interactivity is native
plotly `hovertemplate`, which is what this builder uses.

- `curves`: list of `{'name', 'color', 'peak_alt', 'peak_time', 'hours_up'}`, already
  computed — this module stays a pure reshaper, as its docstring promises.
- One `go.Scatter` per curve, `x=dates`, `y=peak_alt`, `mode='lines+markers'`,
  `connectgaps=False` (NaN months must break the line), `legendgroup=name` so it
  matches the main panel's legend grouping, `line=dict(color=_track_color(...))`
  reusing the existing colour helper (htmlplot.py:46) so a target is the same colour
  in all three panels.
- `customdata` stacked `[peak_time, hours_up]`; hovertemplate along the lines of
  `'%{y:.1f}&deg; at %{customdata[0]} &middot; %{customdata[1]:.1f} h above
  minalt'`, with `hovermode='x unified'` and the same monospace `hoverlabel` as
  `build_figure()` (htmlplot.py:238) so the two panels read alike.
- `fig.add_hline(y=minalt, ...)` when `minalt > 1`, matching htmlplot.py:129.
- `fig.add_vline(x=dates[0], ...)` dashed, annotated "tonight", so the user can see
  where the currently-plotted night sits on the seasonal curve.
- Axes: x `tickformat='%b %Y'`, `dtick='M1'`; y `title='Peak altitude during
  astronomical night [deg]'`, `range=[0, 90]`. `uirevision='skywalker-year'`,
  `height=360`, `autosize`.
- Acceptance: rendering with a hand-built two-curve fixture containing an NaN month
  produces a visible gap and no dropped trace.

### 3. `src/skywalker/webapp.py` — caching and wiring

**`TrackRegistry`** gains `self.year_curves = {}` (name -> the dict from
`year_max_altitudes`), guarded by the existing `self.lock`, plus:

```
def year_series(self, names=None):
    """Cached per-target year curves, computed on first use. Moon excluded."""
```

It iterates `ordered_tracks(names)`, skips `is_moon` and any track whose `coords` is
None (the Moon has neither), collects **all** cache misses into a single
`SkyCoord`, computes them in one `self.walker.year_max_altitudes(...)` call,
memoises the results by name, and returns the list of curve
dicts plus `self.walker.year_dates`. `remove()` should also pop from `year_curves`
so a removed-and-re-added target is recomputed rather than going stale — though with
colours already memoised by name (webapp.py:130) the cheaper and more consistent
choice is to **leave the cache keyed by name and never prune it**, exactly as
`self.colors` is not pruned; the curve for a given name is a pure function of its
coordinates and the site, so it can never go stale. Do that, and note it in the
docstring.

**Layout** (`_serve_layout()`, webapp.py:335): immediately after the existing
`dcc.Loading`/`sw-graph` block and before `sw-controls`, insert

- `html.Button('Show year view', id='sw-year-btn', n_clicks=0)`
- `dcc.Store(id='sw-year-on', data=False)` — a bare bool is safe in a Store; the
  module docstring's warning is about numpy/datetime payloads, which this is not
- `dcc.Loading(children=[dcc.Graph(id='sw-year', figure=go.Figure(), ...)],
  type='default', delay_show=300)` wrapped in `html.Div(id='sw-year-wrap',
  style={'display': 'none'})`

**Callbacks** — both are *consumers* under the module's acyclic rule (webapp.py:9-18):
neither writes anything a mutator reads.

- `_toggle_year`: `Input('sw-year-btn', 'n_clicks')`, `State('sw-year-on', 'data')`
  -> `Output('sw-year-on', 'data')`, `Output('sw-year-wrap', 'style')`,
  `Output('sw-year-btn', 'children')` (flips between "Show year view" / "Hide year
  view"). `prevent_initial_call=True`.
- `_view_year`: `Input('sw-year-on', 'data')`, `Input('sw-table', 'data')`,
  `Input('sw-table', 'selected_rows')`, `Input('sw-table', 'active_cell')` ->
  `Output('sw-year', 'figure')`. Returns `no_update` while the store is False, so
  **nothing is computed until the toggle is on** — this is the whole lazy contract.
  Once on, it derives `_names` and `_focus` exactly as `_view` does (webapp.py:488),
  calls `registry.year_series(_names)`, `htmlplot.build_year_figure(...)`, and
  reuses the existing `_apply_focus(fig, focus)` (webapp.py:265) unchanged — it only
  touches `trace.name`, `trace.opacity` and `line.width`, all of which the new traces
  have.
- Acceptance: with the panel hidden, adding a target triggers no year computation
  (verifiable by a `logger.debug` in `set_year_frames`); with it shown,
  checking/unchecking rows redraws instantly from cache, and clicking a table cell
  dims the other curves just as it does in the main panel.

### 4. `README.md` and `TODO.md`

One short subsection under "Interactive web UI" (README.md:159) describing the year
view and its definition of the y value, and a note that it is computed on demand.

## Performance

The whole feature is **two vectorised astropy calls**, regardless of how many targets
or months are involved: one `get_body('sun', ...)` over the flat 13 x 145 = 1885-epoch
grid (a few hundred ms, once per session), and one broadcast `transform_to()` of an
`(N, 1)` `SkyCoord` against those same 1885 epochs. For N = 30 that is ~57,000
altitude evaluations in a single ERFA call — well under a second. There is no Python
loop over dates or targets anywhere; the reshape is what keeps this cheap.

First toggle therefore costs roughly a second; every redraw afterwards (check,
uncheck, highlight) is a pure dict lookup. If first-toggle latency ever becomes
annoying, the lever is dropping `samples` from 145 to 97 (15-minute cadence), costing
at most ~0.2 deg of peak accuracy.

Explicitly rejected: astroplan's `twilight_evening_astronomical`-style solvers. They
are iterative root-finders, far more expensive per call than one vectorised grid, and
this feature needs a monthly maximum, not an exact twilight timestamp.

## Risks and things to watch

- **`utcoffset` is fixed at the start date** (cli.py:304), so nights half a year away
  are centred on local midnight off by up to an hour where DST applies. A max over a
  24-hour window is insensitive to this; it only matters if the definition later
  changes to something windowed. Worth a comment at the top of `set_year_frames()`.
- **Definition mismatch with the table.** `plotdata.track_summary()` searches the
  peak under `sun_alt < 0` (webapp.py:146) while the year panel uses `< -18`. Month 0
  of the curve will therefore sometimes sit slightly below the table's `Peak alt` for
  the same night. This is correct and intended — the request was explicitly
  astronomical night — but the axis title must say so, and the README note should
  call it out so it does not read as a bug.
- **Polar and high-latitude sites** have months with no astronomical night at all;
  `year_night_mask` is then all-False for that month and every target is NaN there.
  That is the right answer, and `connectgaps=False` renders it honestly.
- **Many targets** make a 13-point line chart busy. The shared `legendgroup` means
  clicking a legend entry in the main panel does *not* hide it in the year panel
  (separate figures), so the table checkboxes remain the real filter. Acceptable;
  revisit only if it bites.
- **A target that is never observable all year** produces an all-NaN curve: it still
  claims a legend entry but draws nothing, which reads as a rendering bug. Cheap fix
  worth doing in the same pass — when every month is NaN, append " (never observable)"
  to the legend name and set the trace to `visible='legendonly'`.
- **The Moon is excluded** — a monthly sample of the Moon's peak altitude is an
  aliasing artefact of the synodic month, not information.

## Verification

1. `pip install -e '.[web]'` in the project venv (the `web` extra already pulls
   `plotly` and `dash`; no new dependency is introduced by this change).
2. `skywalker -f examples/example_file.csv --sitefile examples/sitefilename_example.csv -ns 2019-08-23 --skychart --web --web-open`
3. Confirm the page loads at the same speed as before and the year panel is hidden.
4. Click **Show year view**: the panel appears with one curve per target over
   2019-08 .. 2020-08, monthly ticks, a dashed "tonight" marker on the first point,
   and a red dashed `minalt` line.
5. Sanity-check the physics against a known season — NGC253 (Dec -25) from OPD
   (lat -22.5) should peak around Oct-Nov and bottom out around Apr-May; EG274
   (Dec -39) should peak around May-Jun. Cross-check one month by re-running the CLI
   with `-ns 2019-11-23` and comparing the `Peak alt` in the table to the November
   point (allowing for the -18 deg vs 0 deg difference noted above).
6. Uncheck a row -> its curve disappears from the year panel immediately, no
   recompute pause. Click a table cell -> the other curves dim.
7. Add a new target through the form while the panel is open -> it appears in the
   year panel after a short pause (its own single transform), and toggling the panel
   off and on again is instant.
8. Toggle off -> panel hides, main plot and table unaffected.
