import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from survey.singlecell.plotting import _subset_plot_data, get_plot_data


def make_adata():
    obs = pd.DataFrame(
        {
            "high_q": [True, True, False],
            "ct1": pd.Categorical(["A", "A", "B"]),
        },
        index=["cell_1", "cell_2", "cell_3"],
    )
    adata = ad.AnnData(np.ones((3, 1)), obs=obs)
    adata.obsm["X_umap"] = np.arange(6).reshape(3, 2)
    adata.uns["meta"] = {
        "ct1": pd.DataFrame(
            {"color": ["#ff0000", "#0000ff"]},
            index=["A", "B"],
        )
    }
    return adata


def test_anndata_setting_removes_unused_categories():
    assert ad.settings.remove_unused_categories is True

    subset = _subset_plot_data(make_adata(), {"high_q": True})

    assert subset.obs["ct1"].cat.categories.tolist() == ["A"]


def test_get_plot_data_excludes_sliced_categories_from_color_dict():
    subset = _subset_plot_data(make_adata(), {"high_q": True})

    _, _, color_dict = get_plot_data(subset, color="ct1", basis="umap")

    assert color_dict == {"A": "#ff0000"}


def test_mudata_subset_uses_modality_local_column():
    md = pytest.importorskip("mudata")
    adata = make_adata()
    mdata = md.MuData({"rna": adata.copy(), "xyz": adata.copy()})

    subset = _subset_plot_data(mdata, {"rna:high_q": True})

    assert subset.n_obs == 2
    assert subset["rna"].obs["ct1"].cat.categories.tolist() == ["A"]
