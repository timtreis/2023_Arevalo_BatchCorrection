"""Generate reference.h5ad for a trained scPoli atlas by loading the model
and running get_latent() on the reference data. Run when nb20 failed after
model.save() before the reference embedding step."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.paths import DATA_OUT, MODEL_OUT, scenario_input_parquet
from src.data_io import load_parquet_as_anndata

import anndata as ad
import numpy as np
from scarches.models.scpoli import scPoli

REFERENCE_SCENARIO = "scenario_5"
ATLAS_NAME = f"scpoli_atlas_{REFERENCE_SCENARIO}_full"
BATCH_KEY = "Metadata_Source"
LABEL_KEY = "Metadata_JCP2022"

atlas_dir = MODEL_OUT / ATLAS_NAME
print(f"atlas_dir: {atlas_dir}")

print("loading reference data...")
ref = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
print(f"ref shape: {ref.shape}")

# Reconstruct the _cell_type column (needed by scPoli model architecture)
import json
from src.target2 import load_manifest
t2_ids = load_manifest(DATA_OUT / "target2_compounds.json")
ref.obs["_cell_type"] = ref.obs[LABEL_KEY].astype(str).where(
    ref.obs[LABEL_KEY].isin(t2_ids), "BACKGROUND"
)

print("loading trained scPoli model...")
model = scPoli.load(str(atlas_dir), adata=ref)
print("model loaded")

model.model.eval()
print("computing reference embeddings...")
ref_emb = model.get_latent(ref, mean=True)
print(f"reference embedding shape: {ref_emb.shape}")

ref_out = ad.AnnData(obs=ref.obs[[BATCH_KEY, LABEL_KEY]])
ref_out.obsm["X_scpoli_mapped"] = ref_emb
out_path = atlas_dir / "reference.h5ad"
ref_out.write_h5ad(out_path)
print(f"saved reference embedding to {out_path}")
