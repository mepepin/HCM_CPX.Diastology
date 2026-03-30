from __future__ import annotations

import json
import re
from pathlib import Path


def _is_runtime_example_cell(source: str) -> bool:
    stripped = source.strip()
    if not stripped:
        return True

    hard_skip_patterns = [
        r"pd\.read_csv\s*\(",
        r"pd\.read_excel\s*\(",
        r"\bcsv_path\s*=\s*Path\(",
    ]
    if any(re.search(pattern, source) for pattern in hard_skip_patterns):
        return True

    definition_markers = ("def ", "class ", "@dataclass")
    if any(marker in source for marker in definition_markers):
        return False

    runtime_patterns = [
        r"\bomarx_fit\s*\(",
        r"\bsave_omarx_diagnostics\s*\(",
        r"\bhelper_feature_importance\s*\(",
        r"\bcfg_[A-Za-z0-9_]+\b",
        r"\bbasis_[A-Za-z0-9_]+\b",
        r"\bdf_model\b",
        r"\bdf_vo2_imp\b",
        r"\bdf_vevco2_imp\b",
    ]
    return any(re.search(pattern, source) for pattern in runtime_patterns)


def _load_notebook_namespace() -> dict:
    notebook_path = Path(__file__).resolve().parents[1] / "omarx_v3.ipynb"
    if not notebook_path.exists():
        raise FileNotFoundError(f"OMARX notebook not found: {notebook_path}")

    notebook = json.loads(notebook_path.read_text())
    code_blocks = []

    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue

        source = "".join(cell.get("source", []))
        if not source.strip():
            continue

        if 'if __name__ == "__main__":' in source:
            source = source.split('if __name__ == "__main__":', 1)[0]

        if _is_runtime_example_cell(source):
            continue

        if source.strip():
            code_blocks.append(source)

    namespace = {"__name__": "omarx_from_notebook"}
    exec("\n\n".join(code_blocks), namespace)
    return namespace


_NS = _load_notebook_namespace()

OMARXConfig = _NS["OMARXConfig"]
omarx_fit = _NS["omarx_fit"]
helper_build_basis_for_knots = _NS["helper_build_basis_for_knots"]
helper_fit_model = _NS["helper_fit_model"]
helper_detect_tracks = _NS["helper_detect_tracks"]
helper_standardize = _NS["helper_standardize"]
helper_build_interactions = _NS.get("helper_build_interactions")
helper_feature_importance = _NS.get("helper_feature_importance")
save_omarx_diagnostics = _NS.get("save_omarx_diagnostics")
sm = _NS["sm"]
np = _NS["np"]
pd = _NS["pd"]

__all__ = [
    "OMARXConfig",
    "omarx_fit",
    "helper_build_basis_for_knots",
    "helper_fit_model",
    "helper_detect_tracks",
    "helper_standardize",
    "helper_build_interactions",
    "helper_feature_importance",
    "save_omarx_diagnostics",
    "sm",
    "np",
    "pd",
]
