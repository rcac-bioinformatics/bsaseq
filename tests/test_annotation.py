"""Tests for snpEff annotation functionality."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from bsaseq.analysis.candidates import (
    AnnotatedCandidate,
    CandidateGene,
    CandidateRegion,
    CandidateVariant,
    annotate_candidates,
    summarize_candidate_genes,
)
from bsaseq.annotation.snpeff import (
    VariantAnnotation,
    VariantImpact,
    check_snpeff_available,
    get_best_annotation,
)
from bsaseq.io.vcf import write_candidates_vcf
from bsaseq.io.writers import (
    write_annotated_candidates_tsv,
    write_candidate_genes_tsv,
)


def make_candidate_variant(
    chrom: str = "chr1",
    pos: int = 1000,
    ref: str = "A",
    alt: str = "G",
    af_high: float = 0.95,
    af_low: float = 0.05,
    delta_af: float = 0.9,
    dp_high: int = 50,
    dp_low: int = 50,
    region_rank: int = 1,
    distance_to_peak: int = 0,
) -> CandidateVariant:
    """Helper to create test candidate variants."""
    return CandidateVariant(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        af_high=af_high,
        af_low=af_low,
        delta_af=delta_af,
        dp_high=dp_high,
        dp_low=dp_low,
        region_rank=region_rank,
        distance_to_peak=distance_to_peak,
    )


def make_annotation(
    chrom: str = "chr1",
    pos: int = 1000,
    ref: str = "A",
    alt: str = "G",
    effect: str = "missense_variant",
    impact: VariantImpact = VariantImpact.MODERATE,
    gene_name: str = "GENE1",
    gene_id: str = "GENE1_ID",
) -> VariantAnnotation:
    """Helper to create test annotations."""
    return VariantAnnotation(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        effect=effect,
        impact=impact,
        gene_name=gene_name,
        gene_id=gene_id,
        feature_type="transcript",
        feature_id="TRANS1",
        biotype="protein_coding",
        hgvs_c="c.123A>G",
        hgvs_p="p.Lys41Arg",
        cdna_pos="123",
        cds_pos="123",
        protein_pos="41",
    )


class TestVariantImpact:
    """Tests for VariantImpact enum."""

    def test_impact_values(self) -> None:
        """Test impact enum values."""
        assert VariantImpact.HIGH.value == "HIGH"
        assert VariantImpact.MODERATE.value == "MODERATE"
        assert VariantImpact.LOW.value == "LOW"
        assert VariantImpact.MODIFIER.value == "MODIFIER"

    def test_impact_rank(self) -> None:
        """Test impact ranking (lower = more severe)."""
        assert VariantImpact.HIGH.rank == 0
        assert VariantImpact.MODERATE.rank == 1
        assert VariantImpact.LOW.rank == 2
        assert VariantImpact.MODIFIER.rank == 3

    def test_impact_sorting(self) -> None:
        """Test that impacts sort correctly by rank."""
        impacts = [
            VariantImpact.LOW,
            VariantImpact.HIGH,
            VariantImpact.MODIFIER,
            VariantImpact.MODERATE,
        ]
        sorted_impacts = sorted(impacts, key=lambda i: i.rank)
        assert sorted_impacts == [
            VariantImpact.HIGH,
            VariantImpact.MODERATE,
            VariantImpact.LOW,
            VariantImpact.MODIFIER,
        ]


class TestVariantAnnotation:
    """Tests for VariantAnnotation dataclass."""

    def test_annotation_creation(self) -> None:
        """Test creating a VariantAnnotation."""
        ann = make_annotation()
        assert ann.chrom == "chr1"
        assert ann.pos == 1000
        assert ann.effect == "missense_variant"
        assert ann.impact == VariantImpact.MODERATE
        assert ann.gene_name == "GENE1"

    def test_is_coding(self) -> None:
        """Test is_coding property."""
        # Coding variant
        ann = make_annotation()
        assert ann.is_coding is True

        # Non-coding (empty cds_pos)
        ann2 = VariantAnnotation(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            effect="intron_variant",
            impact=VariantImpact.MODIFIER,
            gene_name="GENE1",
            gene_id="GENE1_ID",
            feature_type="transcript",
            feature_id="TRANS1",
            biotype="protein_coding",
            hgvs_c="",
            hgvs_p="",
            cdna_pos="",
            cds_pos="",
            protein_pos="",
        )
        assert ann2.is_coding is False

    def test_is_lof(self) -> None:
        """Test is_lof property for loss-of-function variants."""
        lof_effects = [
            "stop_gained",
            "stop_lost",
            "start_lost",
            "frameshift_variant",
            "splice_acceptor_variant",
            "splice_donor_variant",
            "transcript_ablation",
        ]

        for effect in lof_effects:
            ann = make_annotation(effect=effect, impact=VariantImpact.HIGH)
            assert ann.is_lof is True, f"{effect} should be LOF"

        # Non-LOF effects
        non_lof = make_annotation(effect="missense_variant")
        assert non_lof.is_lof is False


class TestCheckSnpeffAvailable:
    """Tests for check_snpeff_available function."""

    def test_snpeff_not_available(self) -> None:
        """Test when snpEff is not in PATH."""
        with patch("shutil.which", return_value=None):
            assert check_snpeff_available() is False

    def test_snpeff_available(self) -> None:
        """Test when snpEff is in PATH."""
        with patch("shutil.which", return_value="/usr/local/bin/snpEff"):
            assert check_snpeff_available() is True


class TestGetBestAnnotation:
    """Tests for get_best_annotation function."""

    def test_empty_list(self) -> None:
        """Test with empty list."""
        assert get_best_annotation([]) is None

    def test_single_annotation(self) -> None:
        """Test with single annotation."""
        ann = make_annotation()
        assert get_best_annotation([ann]) == ann

    def test_prefer_high_impact(self) -> None:
        """Test that HIGH impact is preferred over others."""
        high = make_annotation(effect="stop_gained", impact=VariantImpact.HIGH)
        moderate = make_annotation(effect="missense_variant", impact=VariantImpact.MODERATE)
        low = make_annotation(effect="synonymous_variant", impact=VariantImpact.LOW)

        best = get_best_annotation([moderate, low, high])
        assert best == high

    def test_prefer_lof(self) -> None:
        """Test that LOF is preferred when prioritize_lof=True."""
        lof = make_annotation(effect="frameshift_variant", impact=VariantImpact.HIGH)
        non_lof = make_annotation(effect="missense_variant", impact=VariantImpact.HIGH)

        best = get_best_annotation([non_lof, lof], prioritize_lof=True)
        assert best == lof

    def test_no_lof_priority(self) -> None:
        """Test without LOF prioritization."""
        lof = make_annotation(effect="frameshift_variant", impact=VariantImpact.HIGH, gene_name="GENE_A")
        non_lof = make_annotation(effect="stop_gained", impact=VariantImpact.HIGH, gene_name="GENE_B")

        # Both are HIGH impact, should sort by gene name
        best = get_best_annotation([lof, non_lof], prioritize_lof=False)
        # Should prefer alphabetically first gene when equal
        assert best.gene_name == "GENE_A"


class TestWriteCandidatesVcf:
    """Tests for write_candidates_vcf function."""

    def test_write_vcf(self, tmp_path: Path) -> None:
        """Test writing candidates to VCF format."""
        candidates = [
            make_candidate_variant(chrom="chr1", pos=1000, ref="A", alt="G"),
            make_candidate_variant(chrom="chr1", pos=2000, ref="C", alt="T"),
            make_candidate_variant(chrom="chr2", pos=500, ref="G", alt="A"),
        ]

        output = tmp_path / "candidates.vcf"
        write_candidates_vcf(candidates, output)

        assert output.exists()
        content = output.read_text()

        # Check header
        assert "##fileformat=VCFv4.2" in content
        assert "##source=bsaseq" in content
        assert "#CHROM\tPOS\tID\tREF\tALT" in content

        # Check variants
        lines = content.strip().split("\n")
        data_lines = [l for l in lines if not l.startswith("#")]
        assert len(data_lines) == 3

        # Check first variant
        first = data_lines[0].split("\t")
        assert first[0] == "chr1"
        assert first[1] == "1000"
        assert first[3] == "A"
        assert first[4] == "G"
        assert "DELTA_AF=" in first[7]


class TestAnnotatedCandidate:
    """Tests for AnnotatedCandidate dataclass."""

    def test_creation_without_annotation(self) -> None:
        """Test creating AnnotatedCandidate without annotation."""
        candidate = make_candidate_variant()
        annotated = AnnotatedCandidate.from_candidate(candidate, annotation=None)

        assert annotated.chrom == "chr1"
        assert annotated.pos == 1000
        assert annotated.effect is None
        assert annotated.gene_name is None
        assert annotated.is_lof is False

    def test_creation_with_annotation(self) -> None:
        """Test creating AnnotatedCandidate with annotation."""
        candidate = make_candidate_variant()
        annotation = make_annotation(effect="stop_gained", impact=VariantImpact.HIGH)

        annotated = AnnotatedCandidate.from_candidate(candidate, annotation)

        assert annotated.chrom == "chr1"
        assert annotated.pos == 1000
        assert annotated.effect == "stop_gained"
        assert annotated.impact == "HIGH"
        assert annotated.gene_name == "GENE1"
        assert annotated.is_lof is True
        assert annotated.hgvs_c == "c.123A>G"
        assert annotated.hgvs_p == "p.Lys41Arg"


class TestCandidateGene:
    """Tests for CandidateGene dataclass."""

    def test_priority_score(self) -> None:
        """Test priority score calculation."""
        gene = CandidateGene(
            gene_name="GENE1",
            gene_id="GENE1_ID",
            n_variants=3,
            n_lof=1,
            n_high_impact=2,
            n_moderate_impact=1,
            best_region_rank=1,
            variants=[],
        )

        # Score = region_rank * 1000 - n_lof * 100 - n_high * 10 - n_moderate
        expected = 1 * 1000 - 1 * 100 - 2 * 10 - 1
        assert gene.priority_score == expected

    def test_lower_region_rank_better(self) -> None:
        """Test that lower region rank gives better (lower) priority score."""
        gene1 = CandidateGene(
            gene_name="GENE1",
            gene_id="GENE1_ID",
            n_variants=1,
            n_lof=0,
            n_high_impact=1,
            n_moderate_impact=0,
            best_region_rank=1,
            variants=[],
        )
        gene2 = CandidateGene(
            gene_name="GENE2",
            gene_id="GENE2_ID",
            n_variants=1,
            n_lof=0,
            n_high_impact=1,
            n_moderate_impact=0,
            best_region_rank=2,
            variants=[],
        )

        assert gene1.priority_score < gene2.priority_score


class TestAnnotateCandidates:
    """Tests for annotate_candidates function."""

    def test_annotate_matching(self) -> None:
        """Test annotation with matching variants."""
        candidates = [
            make_candidate_variant(chrom="chr1", pos=1000, ref="A", alt="G"),
            make_candidate_variant(chrom="chr1", pos=2000, ref="C", alt="T"),
        ]

        annotations = {
            ("chr1", 1000, "A", "G"): make_annotation(pos=1000),
        }

        result = annotate_candidates(candidates, annotations)

        assert len(result) == 2
        # First should be annotated
        assert result[0].gene_name == "GENE1"
        assert result[0].effect == "missense_variant"
        # Second should not be annotated
        assert result[1].gene_name is None
        assert result[1].effect is None

    def test_annotate_empty_candidates(self) -> None:
        """Test with empty candidate list."""
        result = annotate_candidates([], {})
        assert result == []


class TestSummarizeCandidateGenes:
    """Tests for summarize_candidate_genes function."""

    def test_summarize_single_gene(self) -> None:
        """Test summarizing variants in a single gene."""
        candidate = make_candidate_variant()
        annotation = make_annotation(
            effect="stop_gained",
            impact=VariantImpact.HIGH,
            gene_name="GENE1",
        )
        annotated = AnnotatedCandidate.from_candidate(candidate, annotation)

        genes = summarize_candidate_genes([annotated])

        assert len(genes) == 1
        assert genes[0].gene_name == "GENE1"
        assert genes[0].n_variants == 1
        assert genes[0].n_lof == 1
        assert genes[0].n_high_impact == 1

    def test_summarize_multiple_genes(self) -> None:
        """Test summarizing variants across multiple genes."""
        annotated = []

        # GENE1: 2 variants, 1 LOF
        for i in range(2):
            candidate = make_candidate_variant(pos=1000 + i)
            effect = "frameshift_variant" if i == 0 else "missense_variant"
            impact = VariantImpact.HIGH
            ann = make_annotation(
                pos=1000 + i,
                effect=effect,
                impact=impact,
                gene_name="GENE1",
            )
            annotated.append(AnnotatedCandidate.from_candidate(candidate, ann))

        # GENE2: 1 variant, no LOF
        candidate = make_candidate_variant(pos=3000)
        ann = make_annotation(
            pos=3000,
            effect="synonymous_variant",
            impact=VariantImpact.LOW,
            gene_name="GENE2",
        )
        annotated.append(AnnotatedCandidate.from_candidate(candidate, ann))

        genes = summarize_candidate_genes(annotated)

        assert len(genes) == 2
        # Should be sorted by priority score (GENE1 has more LOF/HIGH)
        gene1 = next(g for g in genes if g.gene_name == "GENE1")
        gene2 = next(g for g in genes if g.gene_name == "GENE2")

        assert gene1.n_variants == 2
        assert gene1.n_lof == 1
        assert gene1.n_high_impact == 2
        assert gene2.n_variants == 1
        assert gene2.n_lof == 0

    def test_summarize_no_annotations(self) -> None:
        """Test with unannotated candidates."""
        candidate = make_candidate_variant()
        annotated = AnnotatedCandidate.from_candidate(candidate, None)

        genes = summarize_candidate_genes([annotated])

        # Should return empty - no genes without gene_name
        assert len(genes) == 0


class TestWriteAnnotatedCandidatesTsv:
    """Tests for write_annotated_candidates_tsv function."""

    def test_write_tsv(self, tmp_path: Path) -> None:
        """Test writing annotated candidates to TSV."""
        candidate = make_candidate_variant()
        annotation = make_annotation()
        annotated = AnnotatedCandidate.from_candidate(candidate, annotation)

        output = tmp_path / "annotated.tsv"
        count = write_annotated_candidates_tsv([annotated], output)

        assert count == 1
        assert output.exists()

        content = output.read_text()
        lines = content.strip().split("\n")

        # Check header
        header = lines[0].split("\t")
        assert "effect" in header
        assert "impact" in header
        assert "gene_name" in header
        assert "is_lof" in header

        # Check data
        data = lines[1].split("\t")
        assert len(data) == len(header)


class TestWriteCandidateGenesTsv:
    """Tests for write_candidate_genes_tsv function."""

    def test_write_tsv(self, tmp_path: Path) -> None:
        """Test writing candidate genes to TSV."""
        gene = CandidateGene(
            gene_name="GENE1",
            gene_id="GENE1_ID",
            n_variants=2,
            n_lof=1,
            n_high_impact=2,
            n_moderate_impact=0,
            best_region_rank=1,
            variants=[],
        )

        output = tmp_path / "genes.tsv"
        count = write_candidate_genes_tsv([gene], output)

        assert count == 1
        assert output.exists()

        content = output.read_text()
        lines = content.strip().split("\n")

        # Check header
        header = lines[0].split("\t")
        assert "gene_name" in header
        assert "n_lof" in header
        assert "best_region_rank" in header

        # Check data
        data = lines[1].split("\t")
        assert data[1] == "GENE1"  # gene_name after rank
