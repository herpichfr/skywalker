# SkyWalker 2.0.0 beta 1

First release since v1.4.5, and the first as an installable package. It is a beta: the web UI
in particular has been exercised by driving its callbacks, far more than by hand in a browser.
Please report problems at https://github.com/herpichfr/skywalker/issues.

> **Disclaimer:** this package was built relying heavily on AI-assisted coding tools. Use it at your own discretion, and check its results before relying on them for observation planning. It is provided without warranty; see LICENSE.

## Breaking changes

- **SkyWalker is now a pip-installable package** with a `skywalker` console command, replacing
  the single `src/skywalker.py` script. Install with `pip install .` (or `bash install.sh
  --install`), then run `skywalker ...` instead of `python src/skywalker.py ...`.
- **`requirements.txt` is gone**; dependencies, including the required `astroplan` fork, are
  declared in `pyproject.toml`.
- **Python 3.9 or newer is required.**
- **Example files moved** into `examples/`.
- **Sexagesimal RA is no longer silently read as degrees.** Hourangle vs. degrees is detected
  from units or range, ambiguous values are rejected, and `--raunit hour|deg` forces a reading.
  Plots of hourangle targets made with older versions should be regenerated (see the README).

## New

- **Interactive web UI (`--web`, needs the `web` extra).** The staralt plot, optional
  skychart and a live target table in the browser: add, remove, select, highlight and edit
  targets in place; load a CSV of targets; switch observatory; move the session to another
  night; set the minimum altitude; copy or download the selection as CSV. Opens in a new
  browser tab by default (`--no-web-open` to disable).
- **Year view** in the web UI: each target's peak altitude or usable hours of astronomical
  night, month by month over a calendar year.
- **Table columns** for peak altitude and time, best airmass, **HoursObs** (hours above the
  minimum altitude in astronomical night) and Moon distance.
- **Light and dark themes** following the browser. The night sky is coloured by the Sun and the
  Moon's illumination, not by the theme; the skychart shows the red below-minimum-altitude ring.
- **Interactive HTML figures (`--savehtml`, needs the `html` extra).**
- **Hover readout** on the matplotlib figure: every target's altitude and airmass at the
  pointer's time, and the Moon's illumination.
- **Civil and nautical twilight** bands on the staralt plot.
- Block durations in the web table accept `h` / `m` / `s` units; a bare number is hours.

## Install

    git clone https://github.com/herpichfr/skywalker.git
    cd skywalker
    python3 -m venv venv && source venv/bin/activate
    pip install '.[web]'

Needs `git` (for the `astroplan` fork) and internet access for name and site lookups. Tested on
Linux only.

Verified with a clean install from a fresh clone via `install.sh --install` on Python 3.13,
which resolved matplotlib 3.11, pandas 3.0, astropy 8.0, plotly 7.1 and dash 4.4.
