"""Programmatic API for building and rendering pedigrees from code.

This module provides a fluent builder API for constructing pedigree diagrams
without needing to write .ped files.

Example::

    from pypedigree.api import PedigreeBuilder

    builder = PedigreeBuilder("FAM1")
    builder.add_individual("father", sex="male", phenotype="affected")
    builder.add_individual("mother", sex="female")
    builder.add_individual("child", sex="male", father="father", mother="mother",
                           phenotypes={"epilepsy": 1, "IQ": 105})
    fig = builder.render(output="family.png")

    # Or convert to a Pedigree object for more control
    pedigree = builder.build()
"""

from __future__ import annotations

import io
from typing import TYPE_CHECKING, Any

from pypedigree.config import PedigreeStyle, resolve_style
from pypedigree.errors import PedigreeError
from pypedigree.models import (
    Adoption,
    Childlessness,
    Family,
    Gender,
    Individual,
    Pedigree,
    Phenotype,
    PhenotypeDisplayMode,
    Pregnancy,
    ReproductiveRole,
    Sex,
    Union,
)
from pypedigree.parser import _strip_display_prefix, build_notation

if TYPE_CHECKING:
    from matplotlib.figure import Figure


_SEX_MAP = {
    "male": Sex.MALE,
    "m": Sex.MALE,
    "1": Sex.MALE,
    "female": Sex.FEMALE,
    "f": Sex.FEMALE,
    "2": Sex.FEMALE,
    "unknown": Sex.UNKNOWN,
    "u": Sex.UNKNOWN,
    "0": Sex.UNKNOWN,
}

_PHENO_MAP = {
    "affected": Phenotype.AFFECTED,
    "2": Phenotype.AFFECTED,
    "unaffected": Phenotype.UNAFFECTED,
    "1": Phenotype.UNAFFECTED,
    "unknown": Phenotype.UNKNOWN,
    "0": Phenotype.UNKNOWN,
}


class PedigreeBuilder:
    """Fluent builder for constructing a Pedigree programmatically.

    Args:
        family_id: Default family ID for added individuals.
    """

    def __init__(self, family_id: str = "FAM1") -> None:
        self._default_family = family_id
        self._individuals: dict[str, dict[str, Any]] = {}
        self._unions: list[dict[str, Any]] = []
        self._phenotype_names: list[str] = []
        self._phenotype_display_modes: dict[str, PhenotypeDisplayMode] = {}
        self._style: PedigreeStyle | None = None

    def add_individual(
        self,
        sample_id: str,
        *,
        sex: str | Sex = "unknown",
        phenotype: str | Phenotype = "unknown",
        father: str | None = None,
        mother: str | None = None,
        family_id: str | None = None,
        deceased: bool = False,
        proband: bool = False,
        carrier: bool = False,
        obligate_carrier: bool = False,
        date_of_birth: str | None = None,
        phenotypes: dict[str, float | int] | None = None,
        annotations: list[str] | None = None,
        pregnancy: Pregnancy | str | None = None,
        adopted: Adoption | str | None = None,
        childless: Childlessness | str | None = None,
        art: ReproductiveRole | str | None = None,
        gender: Gender | str | None = None,
        consultand: bool = False,
    ) -> PedigreeBuilder:
        """Add an individual to the pedigree.

        Args:
            sample_id: Unique identifier for this individual.
            sex: "male", "female", or "unknown" (or Sex enum).
            phenotype: "affected", "unaffected", or "unknown" (or Phenotype enum).
            father: Sample ID of the father (None if unknown).
            mother: Sample ID of the mother (None if unknown).
            family_id: Override the default family ID.
            deceased: Whether this individual is deceased.
            proband: Whether this is the proband.
            carrier: Whether this individual is a carrier. Drawn as a patterned
                section (NSGC 2022 notation, the default) or a half-filled
                symbol (``carrier_notation="legacy"``).
            obligate_carrier: Whether this is an obligate carrier: a sparser
                pattern, or a central dot with the legacy notation.
            date_of_birth: Date of birth string.
            phenotypes: Extra phenotypes as name→value dict.
            annotations: Text lines to display below the individual's ID label.
            pregnancy: "P" (ongoing pregnancy), "SAB" (miscarriage), "TOP"
                (termination), "ECT" (ectopic) or "SB" (stillbirth), or a
                ``Pregnancy`` value.
            adopted: "in" (``father``/``mother`` are the adoptive parents)
                or "out" (they are the biological parents).
            childless: "choice" (no children by choice or for an unknown
                reason) or "infertility".
            art: "D" (gamete donor) or "S" (surrogate).
            gender: Gender identity, "man", "woman" or "nonbinary"; the
                symbol follows it, and ``sex`` is taken as the sex assigned
                at birth (AMAB/AFAB written below when they differ).
            consultand: Point an arrow without "P" at the symbol.

        Returns:
            Self for method chaining.

        Raises:
            PedigreeError: If the individual exists already or a notation
                value is not recognised.
        """
        if sample_id in self._individuals:
            raise PedigreeError(f"Individual '{sample_id}' already added.")

        parsed_sex = sex if isinstance(sex, Sex) else _SEX_MAP.get(str(sex).lower(), Sex.UNKNOWN)
        parsed_pheno = (
            phenotype
            if isinstance(phenotype, Phenotype)
            else _PHENO_MAP.get(str(phenotype).lower(), Phenotype.UNKNOWN)
        )

        self._individuals[sample_id] = {
            "sample_id": sample_id,
            "family_id": family_id or self._default_family,
            "sex": parsed_sex,
            "phenotype": parsed_pheno,
            "father": father,
            "mother": mother,
            "deceased": deceased,
            "proband": proband,
            "carrier": carrier,
            "obligate_carrier": obligate_carrier,
            "date_of_birth": date_of_birth,
            "phenotypes": phenotypes or {},
            "annotations": annotations or [],
            "notation": build_notation(
                pregnancy=pregnancy,
                adopted=adopted,
                childless=childless,
                art=art,
                gender=gender,
                consultand=consultand,
            ),
        }

        # Track phenotype names (strip prefix if present)
        if phenotypes:
            stripped_phenos: dict[str, float | int] = {}
            for name, val in phenotypes.items():
                stripped_name, mode = _strip_display_prefix(name)
                stripped_phenos[stripped_name] = val
                if stripped_name not in self._phenotype_names:
                    self._phenotype_names.append(stripped_name)
                if mode != PhenotypeDisplayMode.AUTO:
                    self._phenotype_display_modes[stripped_name] = mode
            self._individuals[sample_id]["phenotypes"] = stripped_phenos

        return self

    def add_union(
        self,
        partner1: str,
        partner2: str,
        *,
        consanguineous: bool = False,
        broken: bool = False,
    ) -> PedigreeBuilder:
        """Explicitly add a union between two individuals.

        Unions are auto-inferred from parent references, so this is only needed
        for childless couples or to set consanguineous/broken flags.

        Args:
            partner1: Sample ID of first partner.
            partner2: Sample ID of second partner.
            consanguineous: Whether this is a consanguineous union.
            broken: Whether this relationship is dissolved.

        Returns:
            Self for method chaining.
        """
        self._unions.append(
            {
                "partner1": partner1,
                "partner2": partner2,
                "consanguineous": consanguineous,
                "broken": broken,
            }
        )
        return self

    def set_style(self, style: PedigreeStyle) -> PedigreeBuilder:
        """Set a custom rendering style.

        Args:
            style: PedigreeStyle instance.

        Returns:
            Self for method chaining.
        """
        self._style = style
        return self

    def build(self) -> Pedigree:
        """Build a Pedigree object from the added individuals.

        Automatically infers unions from parent references and creates
        phantom parents for individuals with only one known parent.

        Returns:
            A fully constructed Pedigree.
        """
        # Group individuals by family
        families: dict[str, Family] = {}
        for data in self._individuals.values():
            fid = data["family_id"]
            if fid not in families:
                families[fid] = Family(family_id=fid)

            # Normalize phenotype values
            extra_pheno: dict[str, float | None] = {}
            extra_pheno_raw: dict[str, float | None] = {}
            if data["phenotypes"]:
                # Collect all values for normalization
                all_values: dict[str, list[float]] = {}
                for ind_data in self._individuals.values():
                    for pname, pval in ind_data.get("phenotypes", {}).items():
                        if pval is not None:
                            all_values.setdefault(pname, []).append(float(pval))

                for pname, pval in data["phenotypes"].items():
                    if pval is not None:
                        raw = float(pval)
                        extra_pheno_raw[pname] = raw
                        vals = all_values.get(pname, [raw])
                        vmin, vmax = min(vals), max(vals)
                        if vmax > vmin:
                            extra_pheno[pname] = (raw - vmin) / (vmax - vmin)
                        else:
                            extra_pheno[pname] = 1.0 if raw != 0 else 0.0
                    else:
                        extra_pheno[pname] = None
                        extra_pheno_raw[pname] = None

            ind = Individual(
                sample_id=data["sample_id"],
                family_id=fid,
                sex=data["sex"],
                phenotype=data["phenotype"],
                father_id=data["father"],
                mother_id=data["mother"],
                date_of_birth=data["date_of_birth"],
                is_deceased=data["deceased"],
                is_proband=data["proband"],
                is_carrier=data["carrier"],
                is_obligate_carrier=data["obligate_carrier"],
                extra_phenotypes=extra_pheno,
                extra_phenotypes_raw=extra_pheno_raw,
                annotations={str(i): v for i, v in enumerate(data.get("annotations", []))},
                notation=data["notation"],
            )
            families[fid].individuals[data["sample_id"]] = ind

        # Infer phantom parents for individuals with one known parent
        for fam in families.values():
            phantoms_to_add: list[Individual] = []
            for ind in list(fam.individuals.values()):
                if ind.father_id and ind.father_id not in fam.individuals:
                    phantom = Individual(
                        sample_id=ind.father_id,
                        family_id=fam.family_id,
                        sex=Sex.MALE,
                        is_phantom=True,
                    )
                    phantoms_to_add.append(phantom)
                if ind.mother_id and ind.mother_id not in fam.individuals:
                    phantom = Individual(
                        sample_id=ind.mother_id,
                        family_id=fam.family_id,
                        sex=Sex.FEMALE,
                        is_phantom=True,
                    )
                    phantoms_to_add.append(phantom)
            for p in phantoms_to_add:
                if p.sample_id not in fam.individuals:
                    fam.individuals[p.sample_id] = p

        # Infer unions from parent references
        for fam in families.values():
            union_map: dict[tuple[str, str], list[str]] = {}
            for ind in fam.individuals.values():
                if ind.father_id and ind.mother_id:
                    key = (ind.father_id, ind.mother_id)
                    union_map.setdefault(key, []).append(ind.sample_id)

            for (father_id, mother_id), children in union_map.items():
                # Check if explicit union was added
                is_consang = False
                is_broken = False
                for u_data in self._unions:
                    pair = {u_data["partner1"], u_data["partner2"]}
                    if pair == {father_id, mother_id}:
                        is_consang = u_data["consanguineous"]
                        is_broken = u_data["broken"]
                        break

                fam.unions.append(
                    Union(
                        partner_ids=(father_id, mother_id),
                        children_ids=children,
                        is_consanguinous=is_consang,
                        is_broken=is_broken,
                    )
                )

            # Add explicit childless unions
            for u_data in self._unions:
                pair = {u_data["partner1"], u_data["partner2"]}
                existing_pairs = {frozenset(u.partner_ids) for u in fam.unions}
                if frozenset(pair) not in existing_pairs:
                    p1, p2 = u_data["partner1"], u_data["partner2"]
                    if p1 in fam.individuals and p2 in fam.individuals:
                        fam.unions.append(
                            Union(
                                partner_ids=(p1, p2),
                                children_ids=[],
                                is_consanguinous=u_data["consanguineous"],
                                is_broken=u_data["broken"],
                            )
                        )

        return Pedigree(
            families=families,
            phenotype_names=self._phenotype_names,
            phenotype_display_modes=self._phenotype_display_modes,
        )

    def render(
        self,
        output: str | None = None,
        *,
        style: PedigreeStyle | None = None,
        show_labels: bool = True,
        show_legend: bool = True,
        circular: bool = False,
    ) -> Figure:
        """Build the pedigree and render it.

        Args:
            output: Output file path (PNG/SVG/PDF). None for in-memory only.
            style: Override the builder's style.
            show_labels: Whether to show individual ID labels.
            show_legend: Whether to show the phenotype legend.
            circular: If True, use circular (radial) layout.

        Returns:
            The matplotlib Figure.
        """
        from pypedigree.renderer import render_pedigree

        pedigree = self.build()
        render_style = resolve_style(style or self._style, pedigree.phenotype_names)

        return render_pedigree(
            pedigree,
            style=render_style,
            output=output,
            show_labels=show_labels,
            show_legend=show_legend,
            circular=circular,
        )

    def render_interactive_svg(
        self,
        output: str | None = None,
        *,
        style: PedigreeStyle | None = None,
        show_labels: bool = True,
        show_legend: bool = True,
    ) -> str:
        """Build the pedigree and render as interactive SVG.

        Args:
            output: Output file path. None for in-memory only.
            style: Override the builder's style.
            show_labels: Whether to show individual ID labels.
            show_legend: Whether to show the phenotype and carrier legend.

        Returns:
            The SVG content as a string.
        """
        from pypedigree.svg_renderer import render_interactive_svg

        pedigree = self.build()
        render_style = resolve_style(style or self._style, pedigree.phenotype_names)

        return render_interactive_svg(
            pedigree,
            style=render_style,
            output=output,
            show_labels=show_labels,
            show_legend=show_legend,
        )

    def _repr_png_(self) -> bytes:
        """Jupyter notebook PNG display."""
        fig = self.render()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        import matplotlib.pyplot as plt

        plt.close(fig)
        return buf.getvalue()
