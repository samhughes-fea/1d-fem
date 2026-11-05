"""Moved to :mod:`onedfem.interfaces.blade.blade_beam_job`; this name stays for old imports."""

import sys

from onedfem.interfaces.blade import blade_beam_job as _moved

sys.modules[__name__] = _moved
