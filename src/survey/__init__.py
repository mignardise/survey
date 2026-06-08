import anndata as ad

ad.settings.remove_unused_categories = True

from . import singlecell
from . import spatial
