from dataclasses import dataclass, field
from typing import List
from load_profile.blade_load_profile_row import BladeLoadProfileRow

@dataclass
class BladeLoadProfileBin:
    rows: List[BladeLoadProfileRow] = field(default_factory=list)

    def get_by_tsr(self, tsr: str) -> BladeLoadProfileRow:
        for row in self.rows:
            if row.tsr == tsr:
                return row
        raise ValueError(f"No TSR result found for: {tsr}")
