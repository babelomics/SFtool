# sftool/variant_confirmation/candidate_vcf.py

from __future__ import annotations

import os
import re
import tempfile
import gzip
import shutil
import subprocess

from pathlib import Path
from typing import Iterable, List, Sequence, Callable, Dict

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



class CandidateVcfNormalizationError(CandidateVcfError):
    """
    Base exception for candidate VCF normalization failures.
    """

    pass


class CandidateVcfReferenceError(CandidateVcfNormalizationError):
    """
    Raised when the reference genome or its FASTA index is invalid.
    """

    pass


class CandidateVcfCommandError(CandidateVcfNormalizationError):
    """
    Raised when a bcftools normalization command fails.
    """

    pass


class CandidateVcfResultError(CandidateVcfNormalizationError):
    """
    Raised when the normalized VCF cannot be mapped back to candidates.
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


class CandidateVcfNormalizer:
    """
    Normalize a diagnostic candidate VCF with bcftools.

    The normalizer:

    - harmonizes candidate chromosome names against the reference FASTA;
    - normalizes sequence alleles;
    - left-aligns indels;
    - splits multiallelic records;
    - sorts and bgzip-compresses the result;
    - creates a tabix index;
    - updates VariantCandidate.normalized_* fields.

    Candidate IDs must be preserved throughout normalization because they are
    used to associate normalized records with the original domain objects.
    """

    DEFAULT_FILENAME = "diagnostic_candidates.normalized.vcf.gz"

    def __init__(
            self,
            bcftools_path: str | Path,
            runner: Callable = subprocess.run,
    ):
        self.bcftools_path = self._validate_executable_path(
            bcftools_path
        )
        self.runner = runner

    def normalize(
            self,
            raw_vcf_path: str | Path,
            candidates: Sequence[VariantCandidate],
            reference_fasta_path: str | Path,
            output_path: str | Path,
    ) -> Path:
        """
        Normalize a raw candidate VCF and update its VariantCandidate objects.

        Parameters
        ----------
        raw_vcf_path
            Uncompressed candidate VCF generated by CandidateVcfWriter.
        candidates
            Candidates represented in the raw VCF.
        reference_fasta_path
            Reference FASTA corresponding to the candidate assembly.
        output_path
            Final bgzip-compressed normalized VCF path.

        Returns
        -------
        Path
            Path to diagnostic_candidates.normalized.vcf.gz.
        """
        raw_vcf_path = self._validate_raw_vcf_path(
            raw_vcf_path
        )

        output_path = self._validate_output_path(
            output_path
        )

        reference_fasta_path = self._validate_reference_fasta(
            reference_fasta_path
        )

        validated_candidates = self._validate_candidates(
            candidates
        )

        reference_contigs = self._read_reference_contigs(
            reference_fasta_path
        )

        chromosome_mapping = self._build_chromosome_mapping(
            candidates=validated_candidates,
            reference_contigs=reference_contigs,
        )

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        candidate_by_id = {
            candidate.candidate_id: candidate
            for candidate in validated_candidates
        }

        work_directory = output_path.parent

        compatible_vcf_path = (
                work_directory
                / ".diagnostic_candidates.reference-compatible.vcf"
        )

        normalized_unsorted_path = (
                work_directory
                / ".diagnostic_candidates.normalized.unsorted.vcf.gz"
        )

        temporary_output_path = (
                work_directory
                / f".{output_path.name}.tmp.vcf.gz"
        )

        temporary_index_path = Path(
            f"{temporary_output_path}.tbi"
        )

        final_index_path = Path(
            f"{output_path}.tbi"
        )

        temporary_paths = [
            compatible_vcf_path,
            normalized_unsorted_path,
            Path(f"{normalized_unsorted_path}.tbi"),
            temporary_output_path,
            temporary_index_path,
        ]

        try:
            self._write_reference_compatible_vcf(
                raw_vcf_path=raw_vcf_path,
                output_path=compatible_vcf_path,
                chromosome_mapping=chromosome_mapping,
            )

            self._run_bcftools_norm(
                input_path=compatible_vcf_path,
                output_path=normalized_unsorted_path,
                reference_fasta_path=reference_fasta_path,
            )

            self._run_bcftools_sort(
                input_path=normalized_unsorted_path,
                output_path=temporary_output_path,
            )

            self._run_bcftools_index(
                vcf_path=temporary_output_path,
            )

            normalized_records = self._read_normalized_records(
                temporary_output_path
            )

            self._validate_normalized_records(
                normalized_records=normalized_records,
                candidate_by_id=candidate_by_id,
            )

            self._update_candidates(
                normalized_records=normalized_records,
                candidate_by_id=candidate_by_id,
            )

            os.replace(
                temporary_output_path,
                output_path,
            )

            os.replace(
                temporary_index_path,
                final_index_path,
            )

            return output_path

        except CandidateVcfNormalizationError:
            raise

        except OSError as exc:
            raise CandidateVcfNormalizationError(
                "Could not finalize normalized candidate VCF "
                f"{output_path}: {exc}"
            ) from exc

        finally:
            for temporary_path in temporary_paths:
                if temporary_path.exists():
                    try:
                        temporary_path.unlink()
                    except OSError:
                        pass

    def normalize_to_directory(
            self,
            raw_vcf_path: str | Path,
            candidates: Sequence[VariantCandidate],
            reference_fasta_path: str | Path,
            output_directory: str | Path,
            filename: str = DEFAULT_FILENAME,
    ) -> Path:
        """
        Normalize candidates using the standard SFtool output filename.
        """
        if not isinstance(output_directory, (str, Path)):
            raise TypeError(
                "output_directory must be a string or pathlib.Path"
            )

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

        if not filename.endswith(".vcf.gz"):
            raise ValueError(
                "Normalized candidate VCF must use the .vcf.gz extension"
            )

        return self.normalize(
            raw_vcf_path=raw_vcf_path,
            candidates=candidates,
            reference_fasta_path=reference_fasta_path,
            output_path=Path(output_directory) / filename,
        )

    @staticmethod
    def _validate_executable_path(
            bcftools_path: str | Path,
    ) -> Path:
        if not isinstance(bcftools_path, (str, Path)):
            raise TypeError(
                "bcftools_path must be a string or pathlib.Path"
            )

        bcftools_path = Path(
            bcftools_path
        )

        if not bcftools_path.is_file():
            raise FileNotFoundError(
                f"bcftools executable not found: {bcftools_path}"
            )

        if not os.access(
                bcftools_path,
                os.X_OK,
        ):
            raise PermissionError(
                f"bcftools is not executable: {bcftools_path}"
            )

        return bcftools_path

    @staticmethod
    def _validate_raw_vcf_path(
            raw_vcf_path: str | Path,
    ) -> Path:
        if not isinstance(raw_vcf_path, (str, Path)):
            raise TypeError(
                "raw_vcf_path must be a string or pathlib.Path"
            )

        raw_vcf_path = Path(
            raw_vcf_path
        )

        if not raw_vcf_path.is_file():
            raise FileNotFoundError(
                f"Raw candidate VCF not found: {raw_vcf_path}"
            )

        if raw_vcf_path.suffix.lower() != ".vcf":
            raise ValueError(
                "Raw candidate VCF must use the .vcf extension"
            )

        return raw_vcf_path

    @staticmethod
    def _validate_output_path(
            output_path: str | Path,
    ) -> Path:
        if not isinstance(output_path, (str, Path)):
            raise TypeError(
                "output_path must be a string or pathlib.Path"
            )

        output_path = Path(
            output_path
        )

        if not output_path.name.endswith(".vcf.gz"):
            raise ValueError(
                "Normalized candidate VCF output must use "
                "the .vcf.gz extension"
            )

        return output_path

    @staticmethod
    def _validate_reference_fasta(
            reference_fasta_path: str | Path,
    ) -> Path:
        if not isinstance(reference_fasta_path, (str, Path)):
            raise TypeError(
                "reference_fasta_path must be a string or pathlib.Path"
            )

        reference_fasta_path = Path(
            reference_fasta_path
        )

        if not reference_fasta_path.is_file():
            raise CandidateVcfReferenceError(
                f"Reference FASTA not found: {reference_fasta_path}"
            )

        fasta_index_path = Path(
            f"{reference_fasta_path}.fai"
        )

        if not fasta_index_path.is_file():
            raise CandidateVcfReferenceError(
                "Reference FASTA index not found: "
                f"{fasta_index_path}"
            )

        return reference_fasta_path

    @staticmethod
    def _validate_candidates(
            candidates: Sequence[VariantCandidate],
    ) -> List[VariantCandidate]:
        if isinstance(candidates, (str, bytes)):
            raise InvalidCandidateCollectionError(
                "candidates must be a sequence of VariantCandidate objects"
            )

        if not isinstance(candidates, Sequence):
            raise InvalidCandidateCollectionError(
                "candidates must be a sequence of VariantCandidate objects"
            )

        candidates = list(
            candidates
        )

        if not candidates:
            raise InvalidCandidateCollectionError(
                "At least one diagnostic variant candidate is required"
            )

        candidate_ids = set()

        for index, candidate in enumerate(
                candidates,
                start=1,
        ):
            if not isinstance(candidate, VariantCandidate):
                raise InvalidCandidateCollectionError(
                    f"Candidate {index} must be a VariantCandidate"
                )

            if candidate.candidate_id in candidate_ids:
                raise InvalidCandidateCollectionError(
                    "Candidate IDs must be unique"
                )

            candidate_ids.add(
                candidate.candidate_id
            )

        assemblies = {
            candidate.assembly
            for candidate in candidates
        }

        if len(assemblies) != 1:
            raise InvalidCandidateCollectionError(
                "All candidates must use the same reference assembly"
            )

        return candidates

    @staticmethod
    def _read_reference_contigs(
            reference_fasta_path: Path,
    ) -> set[str]:
        fasta_index_path = Path(
            f"{reference_fasta_path}.fai"
        )

        contigs = set()

        try:
            with fasta_index_path.open(
                    "r",
                    encoding="utf-8",
            ) as fasta_index:
                for line_number, line in enumerate(
                        fasta_index,
                        start=1,
                ):
                    line = line.rstrip(
                        "\n"
                    )

                    if not line:
                        continue

                    fields = line.split(
                        "\t"
                    )

                    if len(fields) < 2:
                        raise CandidateVcfReferenceError(
                            "Malformed FASTA index line "
                            f"{line_number}: {line!r}"
                        )

                    contigs.add(
                        fields[0]
                    )

        except OSError as exc:
            raise CandidateVcfReferenceError(
                f"Could not read FASTA index {fasta_index_path}: {exc}"
            ) from exc

        if not contigs:
            raise CandidateVcfReferenceError(
                f"Reference FASTA index is empty: {fasta_index_path}"
            )

        return contigs

    def _build_chromosome_mapping(
            self,
            candidates: Sequence[VariantCandidate],
            reference_contigs: set[str],
    ) -> Dict[str, str]:
        """
        Map candidate chromosome names to contigs in the configured FASTA.

        Exact matches are preferred. Only unambiguous chr-prefix and
        mitochondrial aliases are accepted.
        """
        mapping = {}

        for candidate in candidates:
            chromosome = candidate.chromosome

            if chromosome in mapping:
                continue

            mapping[chromosome] = self._resolve_reference_contig(
                chromosome=chromosome,
                reference_contigs=reference_contigs,
            )

        return mapping

    @staticmethod
    def _resolve_reference_contig(
            chromosome: str,
            reference_contigs: set[str],
    ) -> str:
        if chromosome in reference_contigs:
            return chromosome

        alternatives = []

        if chromosome.startswith("chr"):
            alternatives.append(
                chromosome.removeprefix("chr")
            )
        else:
            alternatives.append(
                f"chr{chromosome}"
            )

        mitochondrial_aliases = {
            "M",
            "MT",
            "chrM",
            "chrMT",
        }

        if chromosome in mitochondrial_aliases:
            alternatives.extend(
                mitochondrial_aliases
            )

        matches = {
            alternative
            for alternative in alternatives
            if alternative in reference_contigs
        }

        if len(matches) == 1:
            return matches.pop()

        if not matches:
            raise CandidateVcfReferenceError(
                "Candidate chromosome "
                f"{chromosome!r} is not present in the reference FASTA"
            )

        raise CandidateVcfReferenceError(
            "Candidate chromosome "
            f"{chromosome!r} maps ambiguously to reference contigs: "
            + ", ".join(sorted(matches))
        )

    @staticmethod
    def _write_reference_compatible_vcf(
            raw_vcf_path: Path,
            output_path: Path,
            chromosome_mapping: Dict[str, str],
    ) -> None:
        """
        Copy the raw VCF while replacing CHROM with the reference-compatible
        contig. Candidate IDs and INFO traceability fields remain unchanged.
        """
        try:
            with raw_vcf_path.open(
                    "r",
                    encoding="utf-8",
            ) as source, output_path.open(
                "w",
                encoding="utf-8",
                newline="\n",
            ) as destination:
                for line_number, line in enumerate(
                        source,
                        start=1,
                ):
                    if line.startswith("#"):
                        if line.startswith("##contig=<ID="):
                            continue

                        if line.startswith("#CHROM"):
                            for contig in dict.fromkeys(
                                    chromosome_mapping.values()
                            ):
                                destination.write(
                                    f"##contig=<ID={contig}>\n"
                                )

                        destination.write(
                            line
                        )
                        continue

                    fields = line.rstrip(
                        "\n"
                    ).split(
                        "\t"
                    )

                    if len(fields) < 8:
                        raise CandidateVcfResultError(
                            "Malformed raw candidate VCF record at line "
                            f"{line_number}"
                        )

                    chromosome = fields[0]

                    if chromosome not in chromosome_mapping:
                        raise CandidateVcfResultError(
                            "Raw VCF contains an unexpected chromosome: "
                            f"{chromosome}"
                        )

                    fields[0] = chromosome_mapping[
                        chromosome
                    ]

                    destination.write(
                        "\t".join(fields) + "\n"
                    )

        except OSError as exc:
            raise CandidateVcfNormalizationError(
                "Could not prepare reference-compatible candidate VCF: "
                f"{exc}"
            ) from exc

    def _run_bcftools_norm(
            self,
            input_path: Path,
            output_path: Path,
            reference_fasta_path: Path,
    ) -> None:
        command = [
            str(self.bcftools_path),
            "norm",
            "--fasta-ref",
            str(reference_fasta_path),
            "--multiallelics",
            "-any",
            "--check-ref",
            "e",
            "--output-type",
            "z",
            "--output",
            str(output_path),
            str(input_path),
        ]

        self._run_command(
            command
        )

    def _run_bcftools_sort(
            self,
            input_path: Path,
            output_path: Path,
    ) -> None:
        command = [
            str(self.bcftools_path),
            "sort",
            "--output-type",
            "z",
            "--output-file",
            str(output_path),
            str(input_path),
        ]

        self._run_command(
            command
        )

    def _run_bcftools_index(
            self,
            vcf_path: Path,
    ) -> None:
        command = [
            str(self.bcftools_path),
            "index",
            "--force",
            "--tbi",
            str(vcf_path),
        ]

        self._run_command(
            command
        )

    def _run_command(
            self,
            command: List[str],
    ) -> None:
        try:
            self.runner(
                command,
                check=True,
                capture_output=True,
                text=True,
            )

        except subprocess.CalledProcessError as exc:
            stderr = (
                exc.stderr.strip()
                if isinstance(exc.stderr, str)
                else ""
            )

            message = (
                    "Candidate VCF normalization command failed: "
                    + " ".join(command)
            )

            if stderr:
                message += f"\n{stderr}"

            raise CandidateVcfCommandError(
                message
            ) from exc

        except OSError as exc:
            raise CandidateVcfCommandError(
                "Could not execute candidate VCF normalization command: "
                + " ".join(command)
            ) from exc

    @staticmethod
    def _read_normalized_records(
            normalized_vcf_path: Path,
    ) -> List[dict]:
        records = []

        try:
            with gzip.open(
                    normalized_vcf_path,
                    "rt",
                    encoding="utf-8",
            ) as normalized_vcf:
                for line_number, line in enumerate(
                        normalized_vcf,
                        start=1,
                ):
                    if line.startswith("#"):
                        continue

                    fields = line.rstrip(
                        "\n"
                    ).split(
                        "\t"
                    )

                    if len(fields) < 8:
                        raise CandidateVcfResultError(
                            "Malformed normalized candidate VCF record "
                            f"at line {line_number}"
                        )

                    candidate_id = fields[2]

                    if not candidate_id or candidate_id == ".":
                        raise CandidateVcfResultError(
                            "Normalized candidate record has no candidate ID "
                            f"at line {line_number}"
                        )

                    records.append(
                        {
                            "candidate_id": candidate_id,
                            "chromosome": fields[0],
                            "position": int(fields[1]),
                            "reference": fields[3],
                            "alternate": fields[4],
                        }
                    )

        except ValueError as exc:
            raise CandidateVcfResultError(
                "Normalized candidate VCF contains an invalid position"
            ) from exc

        except OSError as exc:
            raise CandidateVcfResultError(
                "Could not read normalized candidate VCF "
                f"{normalized_vcf_path}: {exc}"
            ) from exc

        if not records:
            raise CandidateVcfResultError(
                "Normalized candidate VCF contains no records"
            )

        return records

    @staticmethod
    def _validate_normalized_records(
            normalized_records: Sequence[dict],
            candidate_by_id: Dict[str, VariantCandidate],
    ) -> None:
        observed_ids = set()

        for record in normalized_records:
            candidate_id = record["candidate_id"]

            if candidate_id not in candidate_by_id:
                raise CandidateVcfResultError(
                    "Normalized VCF contains an unknown candidate ID: "
                    f"{candidate_id}"
                )

            if candidate_id in observed_ids:
                raise CandidateVcfResultError(
                    "Normalization produced multiple records for candidate "
                    f"{candidate_id}"
                )

            observed_ids.add(
                candidate_id
            )

        missing_ids = (
                set(candidate_by_id)
                - observed_ids
        )

        if missing_ids:
            raise CandidateVcfResultError(
                "Normalized VCF is missing candidates: "
                + ", ".join(sorted(missing_ids))
            )

    @staticmethod
    def _update_candidates(
            normalized_records: Sequence[dict],
            candidate_by_id: Dict[str, VariantCandidate],
    ) -> None:
        for record in normalized_records:
            candidate = candidate_by_id[
                record["candidate_id"]
            ]

            candidate.set_normalized_coordinates(
                chromosome=record["chromosome"],
                position=record["position"],
                reference=record["reference"],
                alternate=record["alternate"],
            )