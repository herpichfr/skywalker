# Handoff — CSV upload, min altitude, HoursObs, and web-UI polish

Session of 2026-09-23 (second part), branch `dev`.
Supersedes the previous handoff, which covered the page title, dark mode and the Moon-coloured
sky (committed as `bc5375e`). Its design notes that still hold are carried forward below.

## What was done

**1. Load CSV.** A `Load CSV` button (`dcc.Upload`, id `sw-csv-upload`) replaces the table's
targets with a file's, keeping the Moon row. Column rules are the CLI's `-f` rules, now shared:
the NAME / BLINIT / BLOCKTIME defaulting and the RA/DEC parsing moved out of
`set_target_list()` into `Skywalker.parse_target_dataframe()`, which both call. RA and DEC are
required; a file missing either, or unreadable, or with no rows, is rejected with the reason in
`sw-status` and the table is untouched. `TrackRegistry.replace_from_dataframe()` builds every
track first and swaps only if at least one is usable. Per-row failures inside a valid file (a
duplicate name, a target that never rises) are skipped and reported, as startup loading does.

**2. Min altitude control.** A `Min alt [deg]` box and `Set` button on the Night row
(`sw-minalt`, `sw-minalt-apply`), a `_mutate` branch, and `TrackRegistry.set_minalt()`. Accepts
0 <= value < 90. It moves the staralt line, the skychart ring, the year view and the HoursObs
column. A site switch keeps it: `apply_site()` / `set_location()` never touch `walker.minalt`.

**3. HoursObs column.** Hours of astronomical night (Sun below -18 deg) with the target above
minalt, computed in `track_summary()` with the same mask and formula as
`year_max_altitudes()`'s `hours_up`. `TrackRegistry.astro_night_mask` is kept alongside
`night_mask` and recomputed in `__init__`, `apply_site()` and `apply_date()`. The Moon row
shows `—`. The column id is `usablehours`; only the header says `HoursObs`.

**4. Table and heading polish.** A visible `SkyWalker - Python observation planner tool`
heading (H2) above `sw-title`. The staralt figure's own plotly title is dropped in the web view
(`_build_figure` sets `title=None`, top margin 30) because it duplicated `sw-title`;
`--savehtml` keeps it. The Moon column header is `MoonDist` (id `moondist`). RA and Dec are
fixed at 95 px.

**5. Year view y axis.** "Hours usable" went wrong after zooming because the y axis shared
the year-only `uirevision`, so a zoom taken under one metric was replayed onto the other's
scale. The y axis now has its own key; see the carried-over note.

**6. Red minalt line missing with `-sc`.** plotly's `add_hline(row=, col=)` skips a subplot
that holds no trace yet, and the line was added before any track. It is now added after the
track loop.

**7. Skychart is black at new Moon.** Its background is now black at alpha
`1 - moon_brightness` over white, as the matplotlib skychart draws it, instead of midnightblue.
The staralt Moon-up band is still midnightblue.

**8. `--web` opens the browser by default.** `--web-open` is now a `BooleanOptionalAction`
defaulting to on, with `--no-web-open` to disable. `run_webapp()` calls
`webbrowser.open_new_tab`, only in the parent process when Werkzeug's reloader is active
(`WERKZEUG_RUN_MAIN`), and maps a wildcard bind (`0.0.0.0`, `::`) to `localhost`.

## Files changed

| File | Change |
|---|---|
| `src/skywalker/cli.py` | `parse_target_dataframe()` factored out of `set_target_list()`; `--web-open` on by default with `--no-web-open` |
| `src/skywalker/plotdata.py` | `track_summary()` takes `astro_mask` / `minalt` and returns the HoursObs value |
| `src/skywalker/htmlplot.py` | minalt line after the track loop; year-view y-axis `uirevision`; black skychart |
| `src/skywalker/webapp.py` | heading; min-alt control; CSV upload; HoursObs / MoonDist columns; RA/Dec widths; no staralt title; `set_minalt()`, `replace_from_dataframe()`, `astro_night_mask`; `open_new_tab` |

## Design notes for whoever picks this up

- **HoursObs and the year view can differ by about 0.1 h.** Same formula, different
  sampling: the table uses the night's own 500-sample track, the year view 145 samples a
  night. The table value is the more accurate one; do not "fix" it to match.
- **`parse_target_dataframe()` is the one place CSV target columns are interpreted.** The
  CLI's `-f` path and the web upload both go through it. The CLI keeps its own RA/DEC presence
  check in `set_target_list()` only to keep its file-named error message.
- **The minalt line depends on a trace existing first.** Anything that moves `add_hline()`
  back above the track loop, or builds a figure with an empty `tracks` list, loses the line
  silently. Today `ordered_tracks()` always includes the Moon, so it is not reachable.
- **The CSV upload is a `_mutate` branch, like every other writer of `sw-table.data`.** Keep
  it that way; see the single-writer note below.

### Carried over from the previous handoffs (still true)

- **`dark=False` must stay the default, and it only restyles the figure chrome.**
  `render()` / `write_html()` / every `--savehtml` path call the builders without it.
  The night bands and skychart background are deliberately *not* theme-dependent in either
  mode — do not reintroduce per-theme band colours; the sky is coloured by the Sun and the Moon.
- **`dcc.Dropdown` in Dash 4 is Dash's own component, themed only through `--Dash-*` custom
  properties** that Dash's bundle sets in one unconditional, light-only `:root` block injected
  after our stylesheet. The dark overrides use `html:root` to win on specificity regardless of
  injection order. Dash's CSS also spells one property two ways (`--Dash-Fill-Inverse-Strong`
  and `--Dash-Fill-Inverse-strong`); both are set on purpose.
- **Inline styles accept `var(--...)`,** including DataTable's `style_*` dicts. That is why the
  table needs no server-side theme branch; keep new colours in the CSS file, not in Python.
- **`sw-theme` is `None` until the clientside callback lands,** and both figure callbacks treat
  anything but `'dark'` as light. The first render on a dark system is therefore light for one
  round trip.
- **The minalt ring is tagged `meta='sw-fixed'`** so `_apply_focus()` never dims it. Any other
  chart furniture added as a trace (rather than a shape) needs the same tag, or highlighting a
  target will dim it.

- **`apply_date()` invalidates the year grid and the curve cache.** `set_year_frames()` derives
  its per-month local midnights from `walker.utcoffset` and `walker.inithour`, both of which
  `set_time()` recomputes from the new date, so the existing grid and every curve on it are
  stale. Same reason `apply_site()` clears them, different input.
- **`apply_date()` must not call `recompute()`**, for the same reason `apply_site()` must not:
  it already holds `TrackRegistry.lock` for the whole rebuild and `threading.Lock` is not
  reentrant. It calls `walker.compute_track()` directly. It also deliberately does not call
  `set_location()` or `set_observer()` — the site is unchanged, only the night.
- **The Dur cell's display must round-trip through its parser.** This is the trap in the
  Dur-cell units change (`bb8be20`).
  Once `30s` is accepted, a formatter with no seconds form renders it `0m`, and the next edit
  of that cell reads `0m` back as zero and silently deletes the block. `track_summary()` now
  formats `%dh%02d` for a whole number of hours and minutes at or above 1 h, `%dm` for whole
  minutes below 1 h, and total seconds with an `s` suffix otherwise — so `2h30m15s` displays as
  `9015s`, which is ugly but reads back exactly. Any change to either side has to preserve that
  property; `_parse_blockdur()` and `track_summary()` are a matched pair and both docstrings
  say so.
- **`track_summary()` rounds the duration to whole seconds before formatting it.** `block_ends`
  comes back from a seconds → astropy-hours → seconds round trip carrying float noise (e.g.
  `5399.999999999996`), which would otherwise pick the wrong format branch. The old code also
  had a latent `int(round(frac * 60))` that could yield 60 and render `1h60`; both are fixed by
  the same rounding.
- **The Add-target form's "Block size [s]" field is still seconds-only** and bypasses
  `_parse_blockdur()` entirely, feeding `_resolve_blocktime()` directly. So `2` in that form is
  2 seconds while `2` in an existing row's Dur cell is 2 hours. This was left alone
  deliberately — the field is labelled `[s]` — but it is now the inconsistent one, and is the
  obvious next thing to unify.
- **`set_year_frames(year=None)` is a no-op once any grid exists.** It does not re-resolve
  `None` to the year of `nightstarts` on every call, only the first. `year_max_altitudes()`
  calls it with no arguments, and under the obvious "always default to nightstarts" reading
  that bare call would rebuild the grid back to the default year and undo the year the user
  just picked in the box.
- **The year-curve cache is a function of the calendar year**, on top of coordinates and site.
  It is keyed by target name alone, so a year switch that does not clear it plots the old
  year's curves against the new axis, silently and plausibly.
- **The year view's `uirevision` is two-tier.** The figure-wide key is the year, so x zoom,
  legend state and trace visibility survive a target or marked-date change and reset on a year
  switch. The y axis has its own key, year *and* metric, so a metric switch always resets it:
  `alt` (fixed 0-90) and `hours` (autorange) have incompatible scales.
- **Exactly one callback writes `sw-table.data` / `selected_rows` / `sw-status` / `sw-title`.**
  The site switch, per-cell edits and now the Recalculate button are all branches of `_mutate`
  rather than additional writers. Do not reach for `allow_duplicate`.
- **The site-exclusivity callback must stay a single callback.** One per direction is a
  circular dependency Dash refuses at registration — the app would not start.
- **`_format_blinit()` accepts exactly `24:00:00`** and rejects everything else out of range;
  `_hours_from_hhmmss()` normalises hour 24 to `0.0` and the past-midnight wrap depends on it.
- **A row's `id` is the target name**, so `active_cell.row_id` is a name. `_resolve_focus()`
  degrades a focus naming no current row to "no focus" rather than dimming the whole figure.

## Verified

Against a live `--web` server on `examples/example_file.csv` (`-sc`), driving
`/_dash-update-component` directly:

- The layout has the H2 heading, a single `Night starts` string, columns ending
  `Best X, HoursObs, MoonDist`, and 95 px RA/Dec widths; the staralt figure has no title.
- HoursObs at minalt 10: EG274 3.8, NGC104 9.1, NGC253 9.1. `Set` to 30: 1.9 / 9.1 / 8.4,
  with the staralt line at 30 and the skychart ring at base 61, r 30. 95 is rejected with the
  table unchanged, and 30 survives a switch to Cerro Paranal.
- A CSV with only NAME and RA returns `bad.csv: DEC column not found.` with the table
  unchanged. A NAME,RA,DEC file replaces the rows, keeps the Moon, and selects all.
- The year view's y `uirevision` is `skywalker-year-<year>-alt` / `-hours`, while the layout
  key stays `skywalker-year-<year>`.
- Sky colour: `rgb(0, 0, 0)` at 0 % illumination, mid grey at 50 %, white at 100 %.
- `parse_args(['--web'])` gives `web_open=True`, `['--web', '--no-web-open']` gives `False`.
  With `webbrowser.open_new_tab` stubbed, a real `main()` run on `--web-host 0.0.0.0` called it
  once with `http://localhost:8783/`.

**Not verified:** no browser check of any of this — the zoom behaviour itself, the new controls
in dark mode, a real browser tab opening. A CSV that includes BLINIT / BLOCKTIME columns was not
uploaded. No automated test suite exists.

## TODO

- If "Hours usable" still misbehaves *without* a metric switch, the likely cause is that an x
  zoom never rescales y to the visible data; that needs a relayout callback.
- Unify the Add-target form's "Block size [s]" field with the Dur cell's unit parsing.
- `TODO.md` still lists "Add individual cell editing to the target table", delivered in
  `81b5c12`. Prune it next time that file is touched.
- Uncommitted, not part of this session's commits: the table column rename `Block` →
  `ObsStart` in `webapp.py`, and the matplotlib hover's Moon `illum. NN%` readout in
  `cli.py` / `hover.py`.
