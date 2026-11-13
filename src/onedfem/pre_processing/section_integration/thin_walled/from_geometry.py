"""Moved to :mod:`onedfem.interfaces.blade.thin_walled_sections`; this name stays for old imports."""

import sys

from onedfem.interfaces.blade import thin_walled_sections as _moved

sys.modules[__name__] = _moved
