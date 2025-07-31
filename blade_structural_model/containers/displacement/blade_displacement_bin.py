from dataclasses import dataclass, field
from typing import List
from displacement.blade_displacement_row import BladeDisplacementRow

@dataclass
class BladeDisplacementBin:
    rows: List[BladeDisplacementRow] = field(default_factory=list)

    def get_by_tsr(self, tsr: str) -> BladeDisplacementRow:
        for row in self.rows:
            if row.tsr == tsr:
                return row
        raise ValueError(f"No TSR result found for: {tsr}")
