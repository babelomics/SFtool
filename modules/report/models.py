# modules/report/models.py

from typing import List, Dict, Any


class ReportTable:
    """
    Logical representation of one report table (e.g. PR, RR, PGx).
    Independent of output format.
    """

    def __init__(self,
                 tab_name: str,
                 rows: List[Dict[str, Any]],
                 metadata: Dict[str, Any] | None = None
                 ):
        self.tab_name = tab_name
        self.rows = rows
        self.metadata = metadata or {}


class SampleReport:
    """
    Container for all report tables belonging to one sample.
    """

    def __init__(self, sample_id: str):
        self.sample_id = sample_id
        self.tables: List[ReportTable] = []

    def add_table(self, table: ReportTable | None):
        if table is None:
            return
        self.tables.append(table)

class CoupleReport:
    """
    RR-only couple-level reproductive risk interpretation.
    """

    def __init__(self,
                 sample_a_id: str,
                 sample_b_id: str,
                 rr_mode: str
                 ):
        self.sample_a_id = sample_a_id
        self.sample_b_id = sample_b_id
        self.rr_mode = rr_mode
        self.tables: List[ReportTable] = []

    def add_table(self, table: ReportTable | None):

        if table is None:
            return
        self.tables.append(table)