# Handoff — night switching and Dur-cell units

Session of 2026-09-21, branch `dev`.
Supersedes the previous handoff, which covered the static Jan–Dec year view (committed as
`255a70f`). Its design notes that still hold are carried forward below.

## What was done

**1. A Recalculate button that moves the session to another night** (committed as `ac1a9d6`).
The `Date` box introduced with the year view only moved that panel's marker line; the staralt
and skychart still showed whatever night the server started with. The box has moved out of
`sw-year-wrap` — hidden until "Show year view" is pressed — into a new always-visible `Night`
row above the Observatory controls, with a `Recalculate` button beside it. Typing still moves
the year marker on its own; the button writes the date into `walker.nightstarts` and rebuilds
the night. The work is in `TrackRegistry.apply_date()`, written as `apply_site()`'s sibling.

**2. The Dur cell now takes units.** A bare number in the table's `Dur` column was read as
seconds, so a user typing `2` for two hours got a 2-second block, displayed as `0m` — reported
as "the number I insert is transformed to a random value always smaller than my input".
`1h`, `1.5h`, `3600s` and `1h30m` were all rejected outright. A bare number is now **hours**,
and `h` / `m` / `s` suffixes are accepted, so `1`, `1h`, `60m` and `3600s` are all the same
3600-second block.

## Files changed

| File | Change |
|---|---|
| `src/skywalker/webapp.py` | `TrackRegistry.apply_date()`; `sw-night-row` layout with `sw-year-date` moved into it and `sw-date-apply` beside it; a `sw-date-apply` branch in `_mutate`; `_parse_blockdur()` rewritten with `_BLOCKDUR_HMS_RE` and `_BLOCKDUR_FORMS` |
| `src/skywalker/plotdata.py` | `track_summary()`'s blockdur formatting gains a seconds form and rounds to whole seconds first |

## Design notes for whoever picks this up

- **`apply_date()` invalidates the year grid and the curve cache.** `set_year_frames()` derives
  its per-month local midnights from `walker.utcoffset` and `walker.inithour`, both of which
  `set_time()` recomputes from the new date, so the existing grid and every curve on it are
  stale. Same reason `apply_site()` clears them, different input.
- **`apply_date()` must not call `recompute()`**, for the same reason `apply_site()` must not:
  it already holds `TrackRegistry.lock` for the whole rebuild and `threading.Lock` is not
  reentrant. It calls `walker.compute_track()` directly. It also deliberately does not call
  `set_location()` or `set_observer()` — the site is unchanged, only the night.
- **The Dur cell's display must round-trip through its parser.** This is the trap in change 2.
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

### Carried over from the previous handoffs (still true)

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

The registered callbacks were pulled out of `app.callback_map`, unwrapped to their closures and
driven directly, with the callback context injected via
`dash._callback_context.context_value.set(AttributeDict(triggered_inputs=[...]))`. There is no
automated test suite in this repo, so none of this is captured as a regression test.

For the Recalculate button, at `lco` from `2024-06-15` to `2024-12-21`:

- The staralt altitude trace genuinely recomputes — M31's peak moves from index 411 to 174 —
  and the skychart's polar `r`/`theta` move with it.
- `utcoffset` −4 h → −3 h, `local_times` now spans the new night, the title reads
  `Night starts: 2024-12-21 @ lco`, and the table's peak column follows (`19.6 @ 07:37` →
  `19.5 @ 20:39`).
- `year_of_frame` is `None` and `year_curves` is empty afterwards.
- Empty, `not-a-date` and `2024-02-30` leave every frame attribute identity-unchanged; a forced
  `set_time()` failure rolls back completely with the tracks intact; a site switch followed by
  a date change neither deadlocks nor leaves a stale curve.

For the Dur cell, by real cell edits through the live callback:

- `1`, `1h`, `60m` and `3600s` all give exactly 3600 s; `2` gives a 2.0 h block displayed
  `2h00`; `1.5` → `1h30`, `0.5` → `30m`, `2h30m15s` → 9015 s.
- Case and whitespace are tolerated (`1H`, `90 M`, ` 1h 30m `).
- Every accepted duration round-trips: the displayed string fed back through a second edit
  yields the identical stored duration, including `30s`, `9015s` and `2h30m15s`.
- `abc`, `-1`, `1x`, `h` and `--` are rejected and leave the block untouched; empty and `—`
  clear it.
- No regression in the neighbouring cells: a `Block` edit keeps the duration, a rename works,
  the Moon row is still protected, and a Recalculate after a Dur edit preserves the block.

**Not verified:** no browser round-trip, for any of this. The `Night` row's CSS and the render
of the new `Date` field and `Recalculate` button are unchecked — run `--web` once.

## TODO

- Unify the Add-target form's "Block size [s]" field with the Dur cell's unit parsing, so `2`
  means the same thing in both. See the design note above.
- `TODO.md` still lists "Add individual cell editing to the target table", which was delivered
  in `81b5c12`. Prune it next time that file is touched.
