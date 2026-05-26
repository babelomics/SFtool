# modules/report/sample_report.py

from typing import List, Dict, Any


class ReportTable:
    """
    Logical representation of one report table (e.g. PR, RR, PGx).
    Independent of output format.
    """

    def __init__(self,
                 name: str,
                 rows: List[Dict[str, Any]],
                 metadata: Dict[str, Any] | None = None
                 ):
        self.tab_name = name
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
