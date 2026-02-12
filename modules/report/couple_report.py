# modules/report/couple_report.py

from typing import List, Dict, Any


class CoupleReport:
    """
    RR-only couple-level reproductive risk interpretation.
    """

    def __init__(self, sample_a_id: str, sample_b_id: str, rr_mode: str):
        self.sample_a_id = sample_a_id
        self.sample_b_id = sample_b_id
        self.rr_mode = rr_mode
        self.tables: List[Dict[str, Any]] = []

    def add_table(self, table):
        self.tables.append(table)
