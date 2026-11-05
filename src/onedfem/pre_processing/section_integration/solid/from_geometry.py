"""Moved to :mod:`onedfem.interfaces.blade.solid_sections`; this name stays for old imports."""

import sys

from onedfem.interfaces.blade import solid_sections as _moved

sys.modules[__name__] = _moved
