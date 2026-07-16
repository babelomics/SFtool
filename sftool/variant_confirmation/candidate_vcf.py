# sftool/variant_confirmation/candidate_vcf.py

from __future__ import annotations

import os
import re
import tempfile

from pathlib import Path
from typing import Iterable, List, Sequence

from sftool.variant_confirmation.models import VariantCandidate


class CandidateVcfError(RuntimeError):
    """
    Base exception for candidate VCF generation failures.
    """

    pass


class InvalidCandidateCollectionError(CandidateVcfError):
    """
    Raised when the candidate collection cannot be written to a VCF.
    """

    pass


class CandidateVcfWriteError(CandidateVcfError):
    """
    Raised when the candidate VCF cannot be written to disk.
    """

    pass


class CandidateVcfWriter:
    """
    Generate an uncompressed VCF containing diagnostic variant candidates.

    The writer uses the original genomic coordinates returned by the
    conversion step. It does not normalize variants, access the reference
    genome, compress the VCF or modify VariantCandidate objects.
    """

    DEFAULT_FILENAME = "diagnostic_candidates.raw.vcf"

    FILE_FORMAT = "VCFv4.2"

    SOURCE = "SFtool"

    ASSEMBLY_HEADER_VALUES = {
        "GRCh37": "GRCh37",
        "GRCh38": "GRCh38",
    }

    _VALID_VCF_ALLELE = re.compile(
        r"^[A-Z]+$"
    )

    def write(
            self,
            candidates: Sequence[VariantCandidate],
            output_path: str | Path,
    ) -> Path:
        """
        Write diagnostic candidates to an uncompressed VCF file.

        Parameters
        ----------
        candidates
            Non-empty sequence of VariantCandidate objects. All candidates
            must belong to the same reference assembly.
        output_path
            Complete path of the VCF file to create.

        Returns
        -------
        Path
            Path of the generated candidate VCF.

        Raises
        ------
        TypeError
            If output_path is not a string or Path.
        InvalidCandidateCollectionError
            If candidates are missing, invalid, duplicated or use different
            reference assemblies.
        CandidateVcfWriteError
            If the output directory or file cannot be written.
        """
        output_path = self._validate_output_path(
            output_path
        )

        validated_candidates = self._validate_candidates(
            candidates
        )

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        header_lines = self._build_header(
            candidates=validated_candidates,
        )

        record_lines = [
            self._build_record(candidate)
            for candidate in validated_candidates
        ]

        self._write_atomically(
            output_path=output_path,
            lines=header_lines + record_lines,
        )

        return output_path

    def write_to_directory(
            self,
            candidates: Sequence[VariantCandidate],
            output_directory: str | Path,
            filename: str = DEFAULT_FILENAME,
    ) -> Path:
        """
        Write the candidate VCF using the standard SFtool filename.

        This convenience method is intended for workflow integration, where
        each sample has its own ``variant_confirmation`` output directory.
        """
        output_directory = self._validate_output_directory(
            output_directory
        )

        filename = self._validate_filename(
            filename
        )

        return self.write(
            candidates=candidates,
            output_path=output_directory / filename,
        )

    def _validate_candidates(
            self,
            candidates: Sequence[VariantCandidate],
    ) -> List[VariantCandidate]:
        """
        Validate the candidate collection before writing any data.
        """
        if isinstance(candidates, (str, bytes)):
            raise InvalidCandidateCollectionError(
                "candidates must be a sequence of VariantCandidate objects"
            )

        if not isinstance(candidates, Sequence):
            raise InvalidCandidateCollectionError(
                "candidates must be a sequence of VariantCandidate objects"
            )

        candidates = list(candidates)

        if not candidates:
            raise InvalidCandidateCollectionError(
                "At least one diagnostic variant candidate is required"
            )

        for index, candidate in enumerate(
                candidates,
                start=1,
        ):
            if not isinstance(candidate, VariantCandidate):
                raise InvalidCandidateCollectionError(
                    "Candidate "
                    f"{index} must be a VariantCandidate, "
                    f"got {type(candidate).__name__}"
                )

            self._validate_candidate_for_vcf(
                candidate=candidate,
                index=index,
            )

        self._validate_single_assembly(
            candidates
        )

        self._validate_unique_candidate_ids(
            candidates
        )

        self._validate_unique_variants(
            candidates
        )

        return candidates

    def _validate_candidate_for_vcf(
            self,
            candidate: VariantCandidate,
            index: int,
    ) -> None:
        """
        Validate fields required to construct one VCF record.

        Most field validation is already performed by VariantCandidate.
        These additional checks enforce requirements specific to VCF output.
        """
        if candidate.assembly not in self.ASSEMBLY_HEADER_VALUES:
            raise InvalidCandidateCollectionError(
                "Candidate "
                f"{index} uses unsupported assembly "
                f"{candidate.assembly!r}"
            )

        if any(
                character in candidate.chromosome
                for character in ("\t", "\n", "\r", " ", ",")
        ):
            raise InvalidCandidateCollectionError(
                "Candidate "
                f"{candidate.candidate_id!r} has an invalid chromosome "
                f"value: {candidate.chromosome!r}"
            )

        if any(
                character in candidate.candidate_id
                for character in ("\t", "\n", "\r", " ", ";")
        ):
            raise InvalidCandidateCollectionError(
                "Candidate ID contains characters that cannot be written "
                f"to the VCF ID field: {candidate.candidate_id!r}"
            )

        self._validate_vcf_allele(
            allele=candidate.reference,
            field_name="reference",
            candidate_id=candidate.candidate_id,
        )

        self._validate_vcf_allele(
            allele=candidate.alternate,
            field_name="alternate",
            candidate_id=candidate.candidate_id,
        )

    def _validate_vcf_allele(
            self,
            allele: str,
            field_name: str,
            candidate_id: str,
    ) -> None:
        """
        Ensure an allele can be represented as a sequence allele in VCF.

        Task 7 supports the SNV/indel sequence alleles returned by the GeneBe
        converter. Symbolic alleles and breakend notation are intentionally
        not accepted.
        """
        if not self._VALID_VCF_ALLELE.fullmatch(
                allele
        ):
            raise InvalidCandidateCollectionError(
                f"Candidate {candidate_id!r} has an invalid "
                f"{field_name} allele for candidate VCF generation: "
                f"{allele!r}"
            )

    @staticmethod
    def _validate_single_assembly(
            candidates: Sequence[VariantCandidate],
    ) -> None:
        """
        Ensure all records in the sample candidate VCF use one assembly.
        """
        assemblies = {
            candidate.assembly
            for candidate in candidates
        }

        if len(assemblies) != 1:
            raise InvalidCandidateCollectionError(
                "All diagnostic variant candidates in one VCF must use "
                "the same reference assembly; received: "
                + ", ".join(sorted(assemblies))
            )

    @staticmethod
    def _validate_unique_candidate_ids(
            candidates: Sequence[VariantCandidate],
    ) -> None:
        """
        Reject duplicate IDs because they break candidate traceability.
        """
        seen = set()
        duplicates = set()

        for candidate in candidates:
            if candidate.candidate_id in seen:
                duplicates.add(
                    candidate.candidate_id
                )

            seen.add(
                candidate.candidate_id
            )

        if duplicates:
            raise InvalidCandidateCollectionError(
                "Candidate IDs must be unique; duplicated IDs: "
                + ", ".join(sorted(duplicates))
            )

    @staticmethod
    def _validate_unique_variants(
            candidates: Sequence[VariantCandidate],
    ) -> None:
        """
        Reject duplicated genomic records.

        GeneBeVariantConverter already removes duplicates, but checking again
        prevents malformed candidate VCFs when the writer is used directly.
        """
        seen = set()
        duplicates = []

        for candidate in candidates:
            variant_key = (
                candidate.chromosome,
                candidate.position,
                candidate.reference,
                candidate.alternate,
                candidate.assembly,
            )

            if variant_key in seen:
                duplicates.append(
                    candidate.get_genomic_variant()
                )

            seen.add(
                variant_key
            )

        if duplicates:
            raise InvalidCandidateCollectionError(
                "Duplicate genomic candidates cannot be written to the VCF: "
                + ", ".join(sorted(set(duplicates)))
            )

    def _build_header(
            self,
            candidates: Sequence[VariantCandidate],
    ) -> List[str]:
        """
        Build a minimal deterministic VCF header.

        Only contigs present in the candidate records are declared. Contig
        lengths are not available at this stage and will not be inferred
        without accessing the reference FASTA.
        """
        assembly = candidates[0].assembly

        contigs = self._get_contigs_in_order(
            candidates
        )

        header = [
            f"##fileformat={self.FILE_FORMAT}",
            f"##source={self.SOURCE}",
            (
                "##reference="
                f"{self.ASSEMBLY_HEADER_VALUES[assembly]}"
            ),
            (
                '##INFO=<ID=ASSEMBLY,Number=1,Type=String,'
                'Description="Reference genome assembly">'
            ),
            (
                '##INFO=<ID=ORIGINAL_VARIANT,Number=1,Type=String,'
                'Description="Original genomic candidate before normalization">'
            ),
        ]

        header.extend(
            f"##contig=<ID={contig}>"
            for contig in contigs
        )

        header.append(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        )

        return header

    @staticmethod
    def _get_contigs_in_order(
            candidates: Sequence[VariantCandidate],
    ) -> List[str]:
        """
        Return unique contigs preserving their first-occurrence order.
        """
        contigs = []
        seen = set()

        for candidate in candidates:
            if candidate.chromosome in seen:
                continue

            seen.add(
                candidate.chromosome
            )

            contigs.append(
                candidate.chromosome
            )

        return contigs

    @staticmethod
    def _build_record(
            candidate: VariantCandidate,
    ) -> str:
        """
        Build one eight-column VCF record.

        The VCF has no sample columns because it represents candidate
        definitions, not observed patient genotypes.
        """
        info = ";".join(
            [
                f"ASSEMBLY={candidate.assembly}",
                (
                    "ORIGINAL_VARIANT="
                    f"{candidate.get_genomic_variant()}"
                ),
            ]
        )

        return "\t".join(
            [
                candidate.chromosome,
                str(candidate.position),
                candidate.candidate_id,
                candidate.reference,
                candidate.alternate,
                ".",
                ".",
                info,
            ]
        )

    @staticmethod
    def _validate_output_path(
            output_path: str | Path,
    ) -> Path:
        """
        Validate and normalize the requested output path.
        """
        if not isinstance(output_path, (str, Path)):
            raise TypeError(
                "output_path must be a string or pathlib.Path"
            )

        output_path = Path(
            output_path
        )

        if not output_path.name:
            raise ValueError(
                "output_path must include a filename"
            )

        if output_path.suffix.lower() != ".vcf":
            raise ValueError(
                "Raw candidate VCF output must use the .vcf extension"
            )

        return output_path

    @staticmethod
    def _validate_output_directory(
            output_directory: str | Path,
    ) -> Path:
        """
        Validate the workflow output directory.
        """
        if not isinstance(output_directory, (str, Path)):
            raise TypeError(
                "output_directory must be a string or pathlib.Path"
            )

        return Path(
            output_directory
        )

    @staticmethod
    def _validate_filename(
            filename: str,
    ) -> str:
        """
        Validate the candidate VCF filename.
        """
        if not isinstance(filename, str):
            raise TypeError(
                "filename must be a string"
            )

        filename = filename.strip()

        if not filename:
            raise ValueError(
                "filename must be a non-empty string"
            )

        if Path(filename).name != filename:
            raise ValueError(
                "filename must not contain directory components"
            )

        if not filename.lower().endswith(".vcf"):
            raise ValueError(
                "Raw candidate VCF filename must use the .vcf extension"
            )

        return filename

    @staticmethod
    def _write_atomically(
            output_path: Path,
            lines: Iterable[str],
    ) -> None:
        """
        Write the complete VCF through a temporary file and atomically replace
        the destination.

        This avoids leaving a partially written VCF when an I/O failure occurs.
        """
        temporary_path = None

        try:
            with tempfile.NamedTemporaryFile(
                    mode="w",
                    encoding="utf-8",
                    newline="\n",
                    prefix=f".{output_path.name}.",
                    suffix=".tmp",
                    dir=output_path.parent,
                    delete=False,
            ) as temporary_file:
                temporary_path = Path(
                    temporary_file.name
                )

                for line in lines:
                    temporary_file.write(
                        line
                    )
                    temporary_file.write(
                        "\n"
                    )

            os.replace(
                temporary_path,
                output_path,
            )

        except OSError as exc:
            if (
                    temporary_path is not None
                    and temporary_path.exists()
            ):
                try:
                    temporary_path.unlink()
                except OSError:
                    pass

            raise CandidateVcfWriteError(
                "Could not write diagnostic candidate VCF "
                f"to {output_path}: {exc}"
            ) from exc