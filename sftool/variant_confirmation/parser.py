# sftool/variant_confirmation/parser.py

from __future__ import annotations

import re

from sftool.variant_confirmation.models import (
    VariantConfirmationRequest,
)


class UnsupportedVariantRepresentationError(ValueError):
    """
    Raised when a diagnostic variant representation cannot be classified.
    """
    pass


class VariantRepresentationParser:
    """
    Detect the representation used for a diagnostic variant.

    Supported representations:

        genomic    chromosome:position:reference:alternate
        hgvsc      accession:c.description
        hgvsg      accession:g.description
        hgvsp      accession:p.description

    This parser performs format classification only. It does not fully
    validate HGVS syntax, accessions, transcripts, coordinates or alleles.
    Full validation and conversion are delegated to GeneBe.
    """

    GENOMIC_PATTERN = re.compile(
        r"""
        ^
        (?P<chromosome>
            (?:chr)?
            (?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)
        )
        :
        (?P<position>[1-9][0-9]*)
        :
        (?P<reference>[ACGTN]+)
        :
        (?P<alternate>[ACGTN]+)
        $
        """,
        re.IGNORECASE | re.VERBOSE,
        )

    HGVSC_PATTERN = re.compile(
        r"""
        ^
        (?P<reference>[^:\s]+)
        :
        c\.
        (?P<description>\S+)
        $
        """,
        re.IGNORECASE | re.VERBOSE,
        )

    HGVSG_PATTERN = re.compile(
        r"""
        ^
        (?P<reference>[^:\s]+)
        :
        g\.
        (?P<description>\S+)
        $
        """,
        re.IGNORECASE | re.VERBOSE,
        )

    HGVSP_PATTERN = re.compile(
        r"""
        ^
        (?P<reference>[^:\s]+)
        :
        p\.
        (?P<description>\S+)
        $
        """,
        re.IGNORECASE | re.VERBOSE,
        )

    def parse(
            self,
            request: VariantConfirmationRequest,
    ) -> VariantConfirmationRequest:
        """
        Detect and assign the representation type to a request.

        Parameters
        ----------
        request
            Variant confirmation request containing the original
            representation.

        Returns
        -------
        VariantConfirmationRequest
            The same request object with representation_type populated.

        Raises
        ------
        TypeError
            If request is not a VariantConfirmationRequest.
        UnsupportedVariantRepresentationError
            If the representation cannot be classified.
        """
        if not isinstance(request, VariantConfirmationRequest):
            raise TypeError(
                "request must be a VariantConfirmationRequest"
            )

        representation_type = self.detect(request.variant)

        request.set_representation_type(
            representation_type
        )

        return request

    def detect(self, variant: str) -> str:
        """
        Return the representation type detected for a variant string.

        HGVS patterns are intentionally permissive. A successful classification
        does not guarantee that the representation is valid or convertible.
        """
        if not isinstance(variant, str):
            raise TypeError(
                "variant must be a string, "
                f"got {type(variant).__name__}"
            )

        variant = variant.strip()

        if not variant:
            raise ValueError(
                "variant must be a non-empty string"
            )

        patterns = (
            ("genomic", self.GENOMIC_PATTERN),
            ("hgvsc", self.HGVSC_PATTERN),
            ("hgvsg", self.HGVSG_PATTERN),
            ("hgvsp", self.HGVSP_PATTERN),
        )

        for representation_type, pattern in patterns:
            if pattern.fullmatch(variant):
                return representation_type

        if variant.count(":") == 3:
            chromosome, position, reference, alternate = variant.split(":")

            if not reference or not alternate:
                raise UnsupportedVariantRepresentationError(
                    "Invalid genomic variant representation: "
                    "REF and ALT must both contain at least one nucleotide. "
                    "Insertions and deletions must include a reference anchor base; "
                    "for example, use '8:342345233:TATC:T' for a deletion or "
                    "'8:342345233:T:TGGA' for an insertion."
                )

        raise UnsupportedVariantRepresentationError(
            "Unsupported diagnostic variant representation: "
            f"{variant!r}. Supported representations are "
            "chr:pos:ref:alt, HGVSc, HGVSg and HGVSp."
)