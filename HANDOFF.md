# Handoff — observatory switcher + editable target table

Session of 2026-09-16, branch `dev`.
Supersedes the previous handoff, which covered the year-view feature (committed as `ad26dde`).

## What was done

Two features added to the `--web` UI, plus five rounds of defect fixes found while verifying them.

**1. Observatory switcher.** A new "Observatory" control row changes the observing site at
runtime, without restarting the server. A searchable dropdown of `EarthLocation.get_site_names()`
or explicit Lat/Lon/Elev boxes with an optional display Name — the two modes are mutually
exclusive in the UI: picking a site clears the coordinate boxes, entering a coordinate clears
the dropdown. Applying a site rebuilds location → observer → time frames → night frames,
clears the year-curve cache, and recomputes every existing target at the new site from its
stored coordinates. Targets and their colours survive the switch; a target that no longer
rises is named in the status line, not silently dropped. Bad input leaves the previous site
fully intact.

**2. Cell-level editing of existing targets.** `Target`, `RA`, `Dec`, `Block` and `Dur` are
editable in place; the computed columns and the Moon row are not. Edits route through the same
resolution and validation path the Add-target form uses. A rename re-keys the track, colour and
year-curve caches; editing a coordinate invalidates that target's curve. An invalid edit reverts
the cell and explains why, so the table never shows a value the registry does not hold. The
Add-target and Remove-selected buttons are unchanged.

## Files changed

| File | Change |
|---|---|
| `src/skywalker/cli.py` | `set_location()` takes explicit `site` / `lat` / `lon` / `elev` / `name`; a bare call is unchanged |
| `src/skywalker/webapp.py` | Observatory row + exclusivity callback; `apply_site()`, `recompute()`, `rename()`, `_resolve_focus()`, `_parse_blockdur()`; per-column `editable`; cell-edit dispatch merged into `_mutate` |
| `src/skywalker/plotdata.py` | `track_blinit_and_blocktime()`; `format_selection_csv()` refactored onto it |

## Design notes for whoever picks this up

- **The year-curve cache is no longer immune to staleness.** The previous handoff and the code
  comments both claimed a curve "cannot go stale" because it is a pure function of coordinates
  and the site. The site can now change, so that is false. `apply_site()` clears the cache
  wholesale; a coordinate edit pops the one entry. The comments have been corrected — do not
  reintroduce the old claim.
- **Exactly one callback writes `sw-table.data` / `selected_rows` / `sw-status` / `sw-title`.**
  Cell editing was merged into the existing `_mutate` rather than added as a second writer, to
  keep the callback graph acyclic. Do not reach for `allow_duplicate`.
- **The site-exclusivity callback must stay a single callback.** One callback per direction
  (dropdown clears coordinates, coordinates clear dropdown) is a circular dependency that Dash
  refuses at registration — the app would not start. It is one callback with all four controls
  as both Input and Output, and it converges because a refire whose *triggering* control is now
  empty returns `no_update` for everything. Removing that guard makes the clearing write bounce
  back and wipe the selection the user just made.
- **`apply_site()` must not call `recompute()`.** It calls `walker.compute_track()` directly.
  `recompute()` takes `TrackRegistry.lock`, `apply_site()` already holds it for the whole
  rebuild, and `threading.Lock` is not reentrant — "DRY-ing" these together deadlocks.
- **`_format_blinit()` accepts exactly `24:00:00`** and rejects everything else out of range.
  `_hours_from_hhmmss()` normalises hour 24 to `0.0` and the past-midnight wrap depends on it,
  so `22:00 → 24:00` and `22:00 → 00:00` both yield a correct 2 h block. The exemption is
  deliberate, not an oversight.
- **A row's `id` is the target name** (`plotdata.track_summary()`), so `active_cell.row_id` is a
  name. A rename, a removal, or a site-switch drop orphans it; `_resolve_focus()` degrades a
  focus naming no current row to "no focus" instead of dimming the whole figure.

## Verified

Both features were exercised by driving the real registered Dash callbacks directly (unwrapped
from `app.callback_map` with an injected callback context): coordinate edits under the default
`raunit='auto'`, curve invalidation, reverts on bad RA / block time / duplicate name, Moon
protection, rename re-keying, site rollback on unknown site and out-of-range latitude and
longitude with `location` identity-unchanged, colour preservation and dropped-target reporting
across a switch, the optional site name and its whitespace fallback, exclusivity convergence in
both directions including `lat=0` at the equator, and stale-focus degradation. `import
skywalker.webapp` and `python -m skywalker --help` pass. There is no automated test suite in
this repo, so none of this is captured as a regression test.

**Not verified:** no browser round-trip. The callback logic is tested but nothing was clicked in
a live page, so CSS and render-level issues in the new Observatory row are unchecked — run
`--web` once. Relatedly, whether `dash_table` bumps `data_timestamp` on a server-driven `data`
write was settled by reading the bundled JS, not by observing a live table; if the "Updated X."
status flashes and vanishes after an edit, that is the thing to revisit.

## TODO

- Add a box to alter the `--nightstarts` date and update the graphs accordingly.
