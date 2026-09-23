# Handoff — page title, dark mode, and Moon-coloured sky

Session of 2026-09-23, branch `dev`.
Supersedes the previous handoff, which covered the Recalculate button and the Dur-cell units
(committed as `ac1a9d6` and `bb8be20`). Its design notes that still hold are carried forward below.

## What was done

**1. Browser tab title.** `Dash(title=...)` is now the fixed string
`SkyWalker - Python observation planner tool` instead of `skywalker <date>`. The night and site
still show in the on-page `sw-title` heading.

**2. `&middot;` rendered literally in the plotly hover panels.** plotly's text renderer decodes
only a handful of HTML entities (`&deg;`, `&amp;`, `&lt;`, `&gt;`, `&nbsp;`, ...); `&middot;` is
not one of them. Every hovertemplate in `htmlplot.py` now uses a literal `·` instead. `&deg;` is
fine and was left alone.

**3. The web UI follows the browser's `prefers-color-scheme`, live.** Two halves:

- *Page chrome by CSS.* New `src/skywalker/assets/skywalker.css` (Dash serves `assets/` next to
  the module automatically; `pyproject.toml` gains `package-data` so non-editable installs ship
  it) defines colour tokens on `:root` with dark overrides under the media query. Every colour
  that used to be an inline literal in `webapp.py` (`_OK_STYLE`, `_ERR_STYLE`, the `_field` label,
  the focus-row highlight, the DataTable header and cells) is now a `var(--...)` reference, so
  Python holds no theme logic for the page and there is no flash on load.
- *Figures by a flag.* Plotly cannot read CSS, so a `dcc.Store(id='sw-theme')` is filled by a
  clientside callback from `matchMedia('(prefers-color-scheme: dark)')`, which also registers a
  `change` listener that pushes updates via `dash_clientside.set_props`. `_view` and `_view_year`
  take it as an Input and pass `dark=` to `htmlplot.build_figure()` / `build_year_figure()`,
  which switch to `plotly_dark` with a `#1e1e1e` paper, a dark hover label, and lighter spike
  and marker lines.

**4. The night sky no longer follows the theme.** First pass of item 3 left the twilight and
night bands translucent, so the theme's plot background showed through them: on a dark page a
bright-Moon night rendered near-black. The bands are now pre-blended over white
(`_over_white()`) and drawn opaque, reproducing the matplotlib figure's colours in both themes.
The Moon-up astronomical night is `midnightblue` at alpha `1 - moon_brightness`, so it is as pale
as the Moon is bright. The polar skychart's whole background takes that same colour, with its
labels and gridlines picking light or dark ink by luminance (`_is_dark()`), and it gains the
matplotlib skychart's red below-`minalt` ring (a full-circle `go.Barpolar`, red at 0.7 over white).

## Files changed

| File | Change |
|---|---|
| `src/skywalker/assets/skywalker.css` | New. Colour tokens, dark overrides, `body`/`input`/`button`, and `--Dash-*` overrides for `dcc.Dropdown` |
| `src/skywalker/webapp.py` | Title; inline colours → `var(--...)`; `sw-theme` Store and clientside callback; `theme` Input on `_view`/`_view_year`; `_build_figure(dark=)`; `_apply_focus` skips `meta='sw-fixed'` traces |
| `src/skywalker/htmlplot.py` | `&middot;` → `·`; `dark=False` on both builders; `_over_white()`, `_rgb_css()`, `_is_dark()`; opaque night bands; Moon-coloured polar background and red minalt ring |
| `pyproject.toml` | `[tool.setuptools.package-data] skywalker = ["assets/*.css"]` |

## Design notes for whoever picks this up

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
- **The skychart and the Moon-up band use `midnightblue`;** the matplotlib skychart uses black.
  At full Moon they look alike; at new Moon the web one is deep blue rather than black. This was
  a choice, not an oversight — flip the one `_over_white('midnightblue', ...)` call in the
  skychart block if a match is wanted.

### Carried over from the previous handoffs (still true)

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
- **`uirevision` is keyed by year, not a constant**, so a metric switch keeps the user's zoom
  and a year switch resets it.
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

Against a live `--web` server on `examples/example_file.csv` (`-sc`), by fetching the layout
and posting to `/_dash-update-component` directly:

- The page `<title>` is the new string; `assets/skywalker.css` is linked and served with 200.
- `sw-theme.data` is registered as a clientside callback on `sw-theme.id`, and is an Input of
  both the `sw-graph` and `sw-year` figure callbacks.
- `theme='light'` returns the default template with no `paper_bgcolor`; `theme='dark'` returns
  `plotly_dark` with `paper_bgcolor='#1e1e1e'` and a `#2a2a2a` hover label.
- The night-band shapes are identical and all opacity 1 in both themes; on that ~94 % Moon night
  the Moon-up band and the polar background are both `rgb(238, 238, 245)`.
- The red ring is `rgb(255, 77, 77)` from `base=91-minalt` over `r=minalt` (10 deg), and keeps full
  opacity while a target is focused.

**Not verified:** no visual check in a browser beyond the user's report that the theme switch
works. That night had no Moon-down stretch, so the opaque black band and the light-ink skychart
branch on a dark-Moon night were not seen live; pick a date near new Moon to check them. No
automated test suite exists, so none of this is a regression test.

## TODO

- Unify the Add-target form's "Block size [s]" field with the Dur cell's unit parsing, so `2`
  means the same thing in both. See the design note above.
- `TODO.md` still lists "Add individual cell editing to the target table", which was delivered
  in `81b5c12`. Prune it next time that file is touched.
- Uncommitted, not part of this session's commit: the table column rename `Block` → `ObsStart`
  in `webapp.py`, and the matplotlib hover's Moon `illum. NN%` readout in `cli.py` / `hover.py`.
