from dataclasses import dataclass, field
from typing import List
from internal_forces.blade_internal_force_row import BladeInternalForceRow

@dataclass
class BladeInternalForceBin:
    rows: List[BladeInternalForceRow] = field(default_factory=list)

    def get_by_tsr(self, tsr: str) -> BladeInternalForceRow:
        for row in self.rows:
            if row.tsr == tsr:
                return row
        raise ValueError(f"No TSR result found for: {tsr}")
