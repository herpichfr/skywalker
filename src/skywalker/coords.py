"""Right Ascension parsing that detects hourangle vs degrees where possible.

RA is the only coordinate with a genuine unit ambiguity: a sexagesimal or bare
decimal value with no explicit units could mean hours (the near-universal
astronomical convention for RA) or degrees. This module resolves that
ambiguity instead of silently guessing, which is what skywalker did before:
``unit=(self.raunit, 'deg')`` with the old ``--raunit`` default of ``'deg'``
read sexagesimal RA as degrees, so ``16:23:33.78`` (16h23m, i.e. 245.89 deg)
was plotted at 16.39 deg -- 89 degrees off the real position.

Declination has no such ambiguity: it is always degrees, so parse_dec() only
has to handle the sexagesimal-vs-decimal formatting, not a unit choice.
"""

import re

from astropy.coordinates import Angle
import astropy.units as u
import numpy as np

_LEADING_NUMBER = re.compile(r'\s*([+-]?\d+(?:\.\d+)?)')


def _leading_number(value):
    """The first numeric field of a sexagesimal or decimal string."""
    match = _LEADING_NUMBER.match(value)
    if match is None:
        raise ValueError(f"Could not parse a number from {value!r}.")
    return float(match.group(1))


def _self_describing_deg(value):
    """Degrees if value carries explicit unit markers (h/d/m/s), else None.

    astropy.coordinates.Angle raises UnitsError when a string has no unit
    markers and none was given, which is exactly the ambiguous case this
    module exists to catch -- so that error is the detector, not a failure.
    """
    try:
        return Angle(value).to(u.deg).value
    except u.UnitsError:
        return None


def _parse_ra_verdict(value, raunit='auto', context=''):
    """Core of parse_ra(): also returns which unit was used ('hour'/'deg').

    Kept separate from parse_ra() so parse_ra_column() can check for a
    column that resolves to more than one unit without parsing every value
    twice.
    """
    _value = str(value).strip()

    if raunit == 'deg':
        return Angle(_value, unit=u.deg).to(u.deg).value, 'deg'
    if raunit == 'hour':
        return Angle(_value, unit=u.hourangle).to(u.deg).value, 'hour'
    if raunit != 'auto':
        raise ValueError(
            f"Invalid raunit {raunit!r}: must be 'auto', 'deg' or 'hour'.")

    try:
        _angle = Angle(_value)               # succeeds only with explicit units
        return _angle.to(u.deg).value, (
            'hour' if _angle.unit == u.hourangle else 'deg')
    except u.UnitsError:
        pass                                  # no unit markers; disambiguate below

    if _leading_number(_value) > 24.:
        return Angle(_value, unit=u.deg).to(u.deg).value, 'deg'

    _as_hour = Angle(_value, unit=u.hourangle).to(u.deg).value
    _as_deg = Angle(_value, unit=u.deg).to(u.deg).value
    _hour_str = Angle(_as_hour, unit=u.deg).to_string(
        unit=u.hourangle, sep='hms', precision=2, pad=True)
    _deg_str = Angle(_as_deg, unit=u.deg).to_string(
        unit=u.deg, sep='dms', precision=2, pad=True)
    raise ValueError(
        f"Ambiguous RA {_value!r}{context}: could be {_hour_str} "
        f"({_as_hour:.3f} deg) or {_deg_str} ({_as_deg:.3f} deg). Pass "
        "--raunit hour or --raunit deg, or write the value with explicit "
        "units (e.g. 16h23m33.78s).")


def parse_ra(value, raunit='auto', context=''):
    """Interpret an RA value, detecting hourangle vs degrees where possible.

    Parameters:
    -----------
    value : str or float
        The RA as given by the user: sexagesimal ("16:23:33.78"), explicit
        ("16h23m33.78s" or "245.89d"), or a bare decimal number.
    raunit : str, optional
        'auto' (default) detects the unit per value; 'hour' or 'deg' force
        every value to that unit, reproducing skywalker's pre-fix behaviour
        for 'deg'.
    context : str, optional
        Appended to the error message, e.g. " for target EG274", to name
        which value was ambiguous.

    Returns the RA in degrees. Raises ValueError when the value is genuinely
    ambiguous and raunit did not say which reading to use.
    """
    _degrees, _verdict = _parse_ra_verdict(value, raunit=raunit,
                                           context=context)
    return _degrees


def parse_ra_column(values, raunit='auto', names=None):
    """Parse a whole RA column, warning if rows resolve to different units.

    Parameters:
    -----------
    values : sequence of str
        One RA value per row.
    raunit : str, optional
        Passed through to parse_ra() for every row. Default is 'auto'.
    names : sequence of str, optional
        Target names, aligned with values, used only to name the row in the
        ambiguity error and in the mixed-unit warning.

    Returns (degrees, warning): degrees is a numpy array; warning is a
    message string when the column resolves to more than one unit (usually a
    sign of a malformed file), or None.
    """
    # list(): a pandas Series indexes by *label*, not position, so a
    # filtered DataFrame (e.g. after --pid) has a non-contiguous index and
    # names[i] would raise KeyError; enumerate() below is always positional.
    values = list(values)
    names = list(names) if names is not None else None

    _degrees = []
    _verdicts = []
    for _i, _value in enumerate(values):
        _name = names[_i] if names is not None else f"row {_i}"
        _deg, _verdict = _parse_ra_verdict(
            _value, raunit=raunit, context=f" for target {_name}")
        _degrees.append(_deg)
        _verdicts.append(_verdict)

    _warning = None
    if len(set(_verdicts)) > 1:
        _warning = ("RA column mixes formats ("
                    + ', '.join(sorted(set(_verdicts)))
                    + "); check the file for a typo.")
    return np.array(_degrees), _warning


def parse_dec(value):
    """Interpret a Dec value in degrees. Always degrees, so no unit choice."""
    return Angle(str(value).strip(), unit=u.deg).to(u.deg).value
