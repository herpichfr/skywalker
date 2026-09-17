# Handoff — static Jan–Dec year view with a date box

Session of 2026-09-17, branch `dev`.
Supersedes the previous handoff, which covered the observatory switcher and in-place target
editing (committed as `81b5c12`). Its design notes that still hold are carried forward below.

## What was done

The `--web` UI's year-view panel was decoupled from `--nightstarts`.

**1. The year grid is a fixed calendar year, January through December.** It was 13
month-spaced nights starting at `--nightstarts`, so the panel's x axis drifted with whatever
night the nightly figure was plotting. It is now 12 nights, the 1st of January through the 1st
of December of one calendar year, defaulting to the year of `--nightstarts`.

**2. A `Date` box in the year panel.** A debounced `YYYY-MM-DD` text input sits next to the
Peak-altitude / Hours-usable radio, defaulting to `--nightstarts`. It is local to the year
view: it never writes `walker.nightstarts`, the nightly figure, or the page title. Choosing a
date in another calendar year re-anchors the whole panel to that year's January–December.

**3. The year figure's vertical line marks the chosen date.** It was a fixed line labelled
`tonight` at `dates[0]`; it is now a dotted line at the box's date, labelled with that date.

## Files changed

| File | Change |
|---|---|
| `src/skywalker/cli.py` | `set_year_frames(year=None, samples=145)` replaces `n_months=13`; builds Jan 1 – Dec 1 of `year`; new `self.year_of_frame` records which year is built |
| `src/skywalker/htmlplot.py` | `build_year_figure()`'s unused `nightstarts` parameter became `marked_date`; the vline moves to it and is labelled with it; `uirevision` keyed by year |
| `src/skywalker/webapp.py` | `_parse_year_date()` + `_YEAR_DATE_RE`; `sw-year-date` input in the panel; `year_series(names, year)` clears the curve cache on a year change; `apply_site()` resets `year_of_frame`; `_view_year` takes the date as a new Input |

## Design notes for whoever picks this up

- **`set_year_frames(year=None)` is a no-op once any grid exists.** It does not resolve `None`
  back to the year of `nightstarts` on every call — only on the first. This is deliberate:
  `year_max_altitudes()` calls `set_year_frames()` with no arguments internally, and under the
  obvious "always default to nightstarts' year" reading that bare call would silently rebuild
  the grid back to the default year and undo the year the user just picked in the box. Only an
  explicit `year` that differs from `self.year_of_frame` triggers a rebuild.
- **The year-curve cache is a function of the calendar year too**, on top of coordinates and
  site. `year_curves` is keyed by target name alone, so a year switch that does not clear it
  plots the old year's curves against the new year's axis — silently, with no error and a
  plausible-looking figure. `year_series()` detects the switch by comparing
  `walker.year_of_frame` before and after `set_year_frames()` and clears the cache wholesale.
  Both the cache's comment and `year_series()`'s docstring have been corrected; do not
  reintroduce the claim that a curve depends only on coordinates and the site.
- **The date box is deliberately not wired to `nightstarts`.** The still-open TODO below asks
  for a box that *does* move the nightly figure. That is a different control; this one was
  specified as independent of it.
- **Malformed or cleared date input returns `no_update`**, keeping whatever the panel already
  shows, rather than falling back to `nightstarts`. A box the user is mid-edit never makes the
  figure jump to a date they did not ask for.
- **`uirevision` is `f'skywalker-year-{dates[0].year}'`, not a constant.** A metric switch or a
  selection change keeps the key, so the user's zoom survives; a year switch changes it and
  resets the zoom, since a range from another year is meaningless against the new axis.
- **The grid's last point is Dec 1**, so a marked date in late December sits past the final
  data point. Plotly autoranges to include the shape, so it should still render — this was not
  confirmed in a browser.

### Carried over from the previous handoff (still true)

- **Exactly one callback writes `sw-table.data` / `selected_rows` / `sw-status` / `sw-title`.**
  Cell editing was merged into `_mutate` rather than added as a second writer, to keep the
  callback graph acyclic. Do not reach for `allow_duplicate`. The year-view date box follows
  the same rule: it is a new `Input` on the existing `_view_year`, not a second writer of
  `sw-year.figure`.
- **The site-exclusivity callback must stay a single callback.** One callback per direction is
  a circular dependency Dash refuses at registration — the app would not start. It is one
  callback with all four controls as both Input and Output, and it converges because a refire
  whose *triggering* control is now empty returns `no_update` for everything.
- **`apply_site()` must not call `recompute()`.** It calls `walker.compute_track()` directly.
  `recompute()` takes `TrackRegistry.lock`, `apply_site()` already holds it for the whole
  rebuild, and `threading.Lock` is not reentrant — "DRY-ing" these together deadlocks.
- **`_format_blinit()` accepts exactly `24:00:00`** and rejects everything else out of range.
  `_hours_from_hhmmss()` normalises hour 24 to `0.0` and the past-midnight wrap depends on it.
  The exemption is deliberate, not an oversight.
- **A row's `id` is the target name** (`plotdata.track_summary()`), so `active_cell.row_id` is
  a name. A rename, a removal, or a site-switch drop orphans it; `_resolve_focus()` degrades a
  focus naming no current row to "no focus" instead of dimming the whole figure.

## Verified

The registered `sw-year.figure` callback was pulled out of `app.callback_map`, unwrapped to its
closure and driven directly — the same technique the previous session used, since there is no
automated test suite in this repo. With `--nightstarts 2024-06-15` at `lco`:

- `walker.year_dates` is 12 dates, `2024-01-01` … `2024-12-01`, `year_of_frame == 2024`.
- A bare `set_year_frames()` after an explicit `set_year_frames(year=2030)` leaves
  `year_of_frame` at 2030 — the no-op semantics above hold.
- The default box value draws the line and its annotation at `2024-06-15`.
- Setting the box to `2031-03-10` rebuilds the grid to Jan–Dec 2031, moves the line, flips
  `uirevision`, and recomputes every curve rather than replaying the 2024 values.
- `''`, `None`, `not-a-date`, `2024-02-30` and `2024-13-01` all return `no_update` with no
  traceback and leave `year_of_frame` alone.
- A date change while the panel is hidden computes nothing.
- `apply_site()` after a year switch completes in ~0.4 s with no deadlock, clears the caches,
  and a subsequent switch to yet another year works.
- Both metrics still render, and `walker.nightstarts` and the page title are untouched
  throughout.
- `import skywalker.webapp` and `python -m skywalker --help` pass.

**Not verified:** no browser round-trip, same gap as the previous session. The new `Date`
field's CSS and render-level appearance are unchecked — run `--web` once.

## TODO

- Add a box to alter the `--nightstarts` date and update the **nightly** graphs accordingly.
  The year view no longer depends on `--nightstarts`, so this is now only about the staralt and
  skychart figures, the page title, and everything derived from `set_time()` / the night frames.
- `TODO.md` still lists "Add individual cell editing to the target table", which was delivered
  in `81b5c12`. Prune it next time that file is touched.
