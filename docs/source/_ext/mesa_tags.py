"""MESA-specific formatting for the generated sphinx-tags overview."""

from __future__ import annotations

from html import escape
from pathlib import Path

import sphinx_tags


ORIGINAL_TAGPAGE = sphinx_tags.tagpage
ORIGINAL_CREATE_FILE = sphinx_tags.Tag.create_file


TAG_GROUPS = (
    (
        "Module",
        (
            "star",
            "binary",
            "astero",
            "adipls",
            "atm",
            "colors",
            "eos",
            "gyre",
            "kap",
            "net",
            "rates",
            "rsp",
            "stella",
        ),
    ),
    (
        "Physics",
        (
            "element-diffusion",
            "gravitational-settling",
            "radiative-levitation",
            "thermohaline",
            "opacity",
            "atmosphere",
            "radiative-gradient",
            "rotation",
            "critical-rotation",
            "gravity-darkening",
            "angular-momentum",
            "conservation",
            "magnetic-field",
            "dynamo",
            "spruit-tayler-dynamo",
            "magnetic-braking",
            "mass-loss",
            "stellar-winds",
            "hot-wind",
            "cool-wind",
            "accretion",
            "eddington",
            "core-helium-burning",
            "hydrogen-burning",
            "helium-burning",
            "carbon-burning",
            "oxygen-burning",
            "silicon-burning",
            "explosive-burning",
            "shell-burning",
            "carbon-flame",
            "conductive-flame",
            "deflagration",
            "flame-speed",
            "carbon-ignition",
            "carbon-flash",
            "off-center-burning",
            "inward-flame",
            "helium-shell-ignition",
            "helium-flash",
            "thermonuclear-runaway",
            "stable-hydrogen-burning",
            "core-hydrogen-burning",
            "triple-alpha",
            "nucleosynthesis",
            "electron-capture",
            "equation-of-state",
            "hydrodynamics",
            "rayleigh-taylor-mixing",
            "starspots",
            "thermal-pulse",
            "c13-pocket",
            "s-process",
            "irradiation",
            "hydrostatic-equilibrium",
            "stability",
        ),
    ),
    (
        "Binaries",
        (
            "binary-evolution",
            "dual-star-evolution",
            "point-mass",
            "mass-transfer",
            "roche-lobe-overflow",
            "explicit-mdot",
            "wind-accretion",
            "black-hole-accretion",
            "gravitational-waves",
            "orbital-angular-momentum",
            "tides",
            "tidal-synchronization",
            "contact-binary",
            "black-hole-binary",
            "high-mass-x-ray-binary",
            "common-envelope",
            "chemically-homogeneous-evolution",
        ),
    ),
    (
        "Numerical methods",
        (
            "composition",
            "entropy",
            "relaxation",
            "brunt-profile",
            "op-mono",
            "type2-opacity",
            "aesopus",
            "t-tau-relation",
            "custom-opacity",
            "mesh-controls",
            "timestep-controls",
            "custom-rates",
            "nuclear-network",
            "adaptive-network",
            "convective-boundary-mixing",
            "convective-premixing",
            "predictive-mixing",
            "partial-mixing-zone",
            "convective-penetration",
            "overshooting",
            "semiconvection",
            "ledoux-criterion",
            "surface-boundary-condition",
            "energy-conservation",
            "split-burn",
            "big-net",
            "burner-substeps",
            "small-timesteps",
            "riemann-solver",
            "hllc",
        ),
    ),
    (
        "Workflows",
        (
            "zams",
            "create-initial-model",
            "prebuilt-model",
            "relax-initial-mass",
            "custom-physics-hooks",
            "save-load",
            "map-output",
            "redo",
            "forced-redo",
            "restart",
            "photo-restart",
            "checksums",
            "null-test",
            "analytic-check",
            "timing",
            "performance",
            "counters",
            "twin-studies",
        ),
    ),
    (
        "Stellar objects and phases",
        (
            "low-mass",
            "intermediate-mass",
            "massive-star",
            "very-massive-star",
            "pre-main-sequence",
            "main-sequence",
            "convective-core",
            "white-dwarf",
            "carbon-oxygen-white-dwarf",
            "helium-white-dwarf",
            "oxygen-neon-white-dwarf",
            "white-dwarf-cooling",
            "compact-remnant",
            "neutron-star",
            "neutron-star-envelope",
            "iron-envelope",
            "hot-subdwarf",
            "red-supergiant",
            "helium-star",
            "wolf-rayet",
            "stellar-envelope",
            "agb",
            "envelope-stripping",
            "ifmr",
            "full-evolution",
            "chandrashekhar-mass",
            "type-ia-progenitor",
            "accretion-induced-collapse",
            "core-collapse",
            "pre-supernova",
            "hydrogen-depletion",
            "supernova",
            "type-iip-supernova",
            "nova",
            "stripped-envelope",
            "pair-instability",
            "pair-instability-supernova",
            "pulsational-pair-instability",
            "low-metallicity",
            "carbon-composition",
            "carbon-oxygen-mixture",
            "kelvin-helmholtz-contraction",
            "black-hole",
            "high-metallicity",
            "zero-metallicity",
            "r-crb-star",
            "hydrogen-deficient",
            "carbon-rich",
            "horizontal-branch",
            "breathing-pulses",
            "very-low-mass",
            "brown-dwarf",
            "brown-dwarf-cooling",
            "planetary-mass-object",
            "irradiated-planet",
            "planetary-cooling",
            "core-envelope-planet",
            "thorne-zytkow-object",
            "solar-model",
            "helium-depletion",
        ),
    ),
    (
        "Pulsation and variability",
        (
            "asteroseismology",
            "p-mode",
            "g-mode",
            "surface-corrections",
            "frequency-matching",
            "seismic-constraints",
            "radial-pulsation",
            "nonradial-pulsation",
            "nonlinear-pulsation",
            "linear-analysis",
            "period-search",
            "instability-strip",
            "cepheid",
            "classical-cepheid",
            "type-ii-cepheid",
            "bl-her",
            "rr-lyrae",
            "delta-scuti",
            "beta-cephei",
            "spb",
            "blap",
            "bep",
            "blue-loop",
            "second-crossing",
            "chaotic-pulsation",
            "helioseismology",
            "pulsation-output",
            "fgong",
        ),
    ),
    (
        "Optimization",
        (
            "optimization",
            "newuoa",
            "simplex",
            "solar-calibration",
            "grid-search",
            "file-input",
            "from-file-search",
            "chi-squared",
            "synthetic-data",
            "custom-parameters",
            "custom-variables",
            "spectroscopic-constraints",
        ),
    ),
)

TAG_GROUP_DISPLAY_ORDER = (
    "Module",
    "Stellar objects and phases",
    "Physics",
    "Pulsation and variability",
    "Binaries",
    "Numerical methods",
    "Workflows",
    "Optimization",
)


def _underline(text: str, character: str) -> str:
    return character * len(text)


def _doc_target(path: str) -> str:
    return str(Path(path).with_suffix(""))


def _tag_count_sort_key(tag: object) -> tuple[int, str, str]:
    return (-len(tag.items), tag.name.lower(), tag.name)


def _group_tags(tags: list[object]) -> list[tuple[str, list[object]]]:
    tags_by_name = {tag.name: tag for tag in tags}
    grouped: list[tuple[str, list[object]]] = []
    used: set[str] = set()
    group_order = {title: index for index, title in enumerate(TAG_GROUP_DISPLAY_ORDER)}

    for title, names in TAG_GROUPS:
        group = sorted(
            [tags_by_name[name] for name in names if name in tags_by_name],
            key=_tag_count_sort_key,
        )
        if group:
            grouped.append((title, group))
            used.update(tag.name for tag in group)

    grouped.sort(key=lambda entry: group_order.get(entry[0], len(group_order)))

    other_tags = sorted(
        [tag for tag in tags if tag.name not in used],
        key=_tag_count_sort_key,
    )
    if other_tags:
        grouped.append(("Other tags", other_tags))

    return grouped


def _format_count(count: int) -> str:
    if count == 1:
        return "1 example"
    return f"{count} examples"


def _tag_card(tag: object) -> str:
    name = escape(tag.name)
    count = len(tag.items)
    count_text = escape(_format_count(count))
    href = escape(f"{tag.file_basename}.html")
    return (
        f'<a class="mesa-tag-card" href="{href}" '
        f'data-tag-name="{name}" data-tag-count="{count}">'
        f'<span class="mesa-tag-name">{name}</span>'
        f'<span class="mesa-tag-count">{count_text}</span>'
        "</a>"
    )


def _raw_html_block(lines: list[str]) -> list[str]:
    return [".. raw:: html", "", *[f"   {line}" for line in lines], ""]


def _overview_controls(tags: list[object]) -> list[str]:
    return _raw_html_block(
        [
            '<div class="mesa-tag-overview" data-mesa-tag-overview>',
            '  <div class="mesa-tag-toolbar" role="group" aria-label="Sort tags">',
            '    <span class="mesa-tag-toolbar-label">Sort tags</span>',
            '    <button class="mesa-tag-sort" type="button" data-mesa-tag-sort="alpha" aria-pressed="false">A-Z</button>',
            '    <button class="mesa-tag-sort is-active" type="button" data-mesa-tag-sort="count" aria-pressed="true">Example count</button>',
            "  </div>",
            "</div>",
        ]
    )


def _overview_script() -> list[str]:
    return _raw_html_block(
        [
            "<script>",
            "(function () {",
            "  function sortTagList(list, mode) {",
            "    var cards = Array.prototype.slice.call(list.querySelectorAll('.mesa-tag-card'));",
            "    cards.sort(function (left, right) {",
            "      var leftName = left.dataset.tagName.toLowerCase();",
            "      var rightName = right.dataset.tagName.toLowerCase();",
            "      if (mode === 'count') {",
            "        var countDiff = Number(right.dataset.tagCount) - Number(left.dataset.tagCount);",
            "        if (countDiff !== 0) return countDiff;",
            "      }",
            "      return leftName.localeCompare(rightName);",
            "    });",
            "    cards.forEach(function (card) { list.appendChild(card); });",
            "  }",
            "",
            "  function setSort(mode) {",
            "    document.querySelectorAll('[data-mesa-tag-list]').forEach(function (list) {",
            "      sortTagList(list, mode);",
            "    });",
            "    document.querySelectorAll('[data-mesa-tag-sort]').forEach(function (button) {",
            "      var isActive = button.dataset.mesaTagSort === mode;",
            "      button.classList.toggle('is-active', isActive);",
            "      button.setAttribute('aria-pressed', isActive ? 'true' : 'false');",
            "    });",
            "  }",
            "",
            "  document.querySelectorAll('[data-mesa-tag-sort]').forEach(function (button) {",
            "    button.addEventListener('click', function () { setSort(button.dataset.mesaTagSort); });",
            "  });",
            "  setSort('count');",
            "}());",
            "</script>",
        ]
    )


def _write_grouped_tagpage(tags, outdir, title, extension, tags_index_head):
    if "rst" not in extension:
        return ORIGINAL_TAGPAGE(tags, outdir, title, extension, tags_index_head)

    all_tags = list(tags.values())
    content = [
        ":orphan:",
        "",
        ".. _tagoverview:",
        "",
        title,
        _underline(title, "#"),
        "",
    ]
    content.extend(_overview_controls(all_tags))

    for group_title, group_tags in _group_tags(all_tags):
        content.extend(
            [
                group_title,
                _underline(group_title, "-"),
                "",
            ]
        )
        content.extend(
            _raw_html_block(
                [
                    '<div class="mesa-tag-grid" data-mesa-tag-list>',
                    *[f"  {_tag_card(tag)}" for tag in group_tags],
                    "</div>",
                ]
            )
        )

    content.extend(_overview_script())

    Path(outdir, "tagsindex.rst").write_text("\n".join(content), encoding="utf8")


def _write_link_tag_file(
    self,
    items,
    extension,
    tags_output_dir,
    srcdir,
    tags_page_title,
    tags_page_header,
):
    if "rst" not in extension:
        return ORIGINAL_CREATE_FILE(
            self,
            items,
            extension,
            tags_output_dir,
            srcdir,
            tags_page_title,
            tags_page_header,
        )

    header = f"{tags_page_title}: {self.name}"
    content = [
        ":orphan:",
        "",
        f".. _sphx_tag_{self.file_basename}:",
        "",
        header,
        _underline(header, "#"),
        "",
        f"{tags_page_header}:",
        "",
    ]

    for item in sorted(items, key=lambda entry: entry.relpath(srcdir).lower()):
        path = item.relpath(srcdir)
        name = Path(path).stem
        content.append(f"- :doc:`{name} <../{_doc_target(path)}>`")

    content.append("")
    Path(srcdir, tags_output_dir, f"{self.file_basename}.rst").write_text(
        "\n".join(content),
        encoding="utf8",
    )


def setup(app):
    sphinx_tags.tagpage = _write_grouped_tagpage
    sphinx_tags.Tag.create_file = _write_link_tag_file

    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
