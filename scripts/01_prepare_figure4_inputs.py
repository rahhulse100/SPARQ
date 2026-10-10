#!/usr/bin/env python3
"""Stage canonical, immutable Figure 4 inputs and audit pathology alignment."""
from __future__ import annotations
import argparse, hashlib, json, shutil, tarfile
from pathlib import Path
import numpy as np
import pandas as pd

ATLAS = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
HER2 = ["A1", "B1", "C1", "D1", "E1", "F1", "G2", "H1"]

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()

def stage_atlas(root: Path, out: Path, manifest: list[dict]):
    source=root/"breast_cancer_atlas"; staged=out/"source_archives"/"atlas"
    staged.mkdir(parents=True,exist_ok=True)
    for name in ["filtered_count_matrices.tar.gz","raw_count_matrices.tar.gz","metadata.tar.gz","spatial.tar.gz"]:
        src=source/name
        if not src.exists(): raise FileNotFoundError(src)
        dst=staged/name
        if not dst.exists(): shutil.copy2(src,dst)
        manifest.append({"cohort":"atlas","source":str(src),"staged":str(dst),"sha256":sha256(dst)})
    for archive in staged.glob("*.tar.gz"):
        with tarfile.open(archive) as tf: tf.extractall(staged)
    for sid in ATLAS:
        base=staged/"filtered_count_matrices"/f"{sid}_filtered_count_matrix"
        spatial=staged/"spatial"/f"{sid}_spatial"
        meta=staged/"metadata"/f"{sid}_metadata.csv"
        needed=[base/"matrix.mtx.gz",base/"barcodes.tsv.gz",base/"features.tsv.gz",spatial/"tissue_positions_list.csv",meta]
        if not all(p.exists() for p in needed): raise RuntimeError(f"Incomplete atlas section {sid}: {needed}")
        sec=out/"canonical"/f"atlas_{sid}"; sec.mkdir(parents=True,exist_ok=True)
        for p in needed: shutil.copy2(p,sec/p.name)
        m=pd.read_csv(meta,index_col=0); m.index=m.index.astype(str)
        pos=pd.read_csv(spatial/"tissue_positions_list.csv",header=None,index_col=0)
        pos.index=pos.index.astype(str); pos.columns=["in_tissue","array_row","array_col","pixel_row","pixel_col"]
        spots=m.join(pos,how="inner"); spots.index.name="spot_id"; spots.reset_index().to_csv(sec/"spots.tsv.gz",sep="\t",index=False,compression="gzip")
        pd.DataFrame({"spot_id":spots.index,"pathology_label":spots["Classification"].astype(str)}).to_csv(sec/"pathology_labels.tsv.gz",sep="\t",index=False,compression="gzip")
        (sec/"section.json").write_text(json.dumps({"cohort":"atlas","section":sid,"n_spots":len(spots)},indent=2))

def stage_her2(root: Path, out: Path, manifest: list[dict]):
    source=root/"figure4_inputs"/"HER2_positive_breast"
    original=root/"HER2_positive_breast"
    archive_dir=out/"source_archives"/"her2"; archive_dir.mkdir(parents=True,exist_ok=True)
    for name in ["count-matrices.zip","images.zip","meta.zip","spot-selections.zip"]:
        src=original/name
        if src.exists():
            dst=archive_dir/name
            if not dst.exists(): shutil.copy2(src,dst)
            manifest.append({"cohort":"her2","section":"all","kind":"source_archive","source":str(src),"staged":str(dst),"sha256":sha256(dst)})
    for sid in HER2:
        count=source/"counts"/"count-matrices"/f"{sid}.tsv.gz"; label=source/"annotations"/f"{sid}_labeled_coordinates.tsv"
        if not (count.exists() and label.exists()): raise RuntimeError(f"Missing HER2 input for {sid}")
        sec=out/"canonical"/f"her2_{sid}"; sec.mkdir(parents=True,exist_ok=True)
        shutil.copy2(count,sec/"counts.tsv.gz"); shutil.copy2(label,sec/"pathology_coordinates.tsv")
        labels=pd.read_csv(label,sep="\t")
        sel_path=source/"selections"/f"{sid}_selection.tsv.gz"
        if not sel_path.exists(): raise RuntimeError(f"Missing HER2 selection file for {sid}: {sel_path}")
        shutil.copy2(sel_path, sec/"selection.tsv.gz")
        selection=pd.read_csv(sel_path,sep="\t",compression="gzip")
        counts=pd.read_csv(count,sep="\t",compression="gzip",usecols=[0],skiprows=1,header=None,names=["count_spot_id"],dtype=str)
        xy=counts.count_spot_id.str.extract(r"^(?P<x>\d+)x(?P<y>\d+)$").astype(int)
        counts=pd.concat([counts,xy],axis=1)
        selection[["x","y"]]=selection[["x","y"]].astype(int)
        direct=counts.merge(selection,left_on=["x","y"],right_on=["x","y"],how="left",validate="one_to_one")
        if direct.new_x.isna().any(): raise RuntimeError(f"{sid}: count-to-selection mapping failed")
        # Iterative one-to-one affine registration: selection new_x/new_y -> pathology x/y.
        # Eight signed coordinate starts prevent a silent x/y reflection mistake.
        from scipy.optimize import linear_sum_assignment
        src=direct[["new_x","new_y"]].to_numpy(float); tgt=labels[["x","y"]].to_numpy(float)
        if len(src)!=len(tgt): raise RuntimeError(f"{sid}: selection/pathology sizes differ")
        def norm(a):
            centered=a-a.mean(0)
            scale=np.sqrt(np.mean((centered**2).sum(1)))
            return centered / (scale if scale > 0 else 1)
        s,t=norm(src),norm(tgt)
        transforms=[]
        for swap in (False, True):
            for sx in (-1, 1):
                for sy in (-1, 1):
                    z=s[:, ::-1] if swap else s.copy()
                    z=z*np.array([sx,sy])
                    transforms.append((f"swap={int(swap)};sx={sx};sy={sy}",z))
        candidates=[]
        for orientation, z in transforms:
            coef=np.linalg.lstsq(np.c_[z,np.ones(len(z))],t,rcond=None)[0]
            for _ in range(100):
                pred=np.c_[z,np.ones(len(z))]@coef
                d=((pred[:,None]-t[None,:])**2).sum(2)**.5
                ii,jj=linear_sum_assignment(d)
                new=np.linalg.lstsq(np.c_[z[ii],np.ones(len(ii))],t[jj],rcond=None)[0]
                if np.max(np.abs(new-coef))<1e-10:
                    coef=new
                    break
                coef=new
            pred=np.c_[z,np.ones(len(z))]@coef
            d=((pred[:,None]-t[None,:])**2).sum(2)**.5
            ii,jj=linear_sum_assignment(d)
            candidates.append((float(d[ii,jj].sum()), orientation, z, coef, ii, jj, d[ii,jj]))
        _, orientation, z, coef, ii, jj, residual=min(candidates,key=lambda x:x[0])
        mapped=direct.iloc[ii].copy().reset_index(drop=True)
        mapped["pathology_spot_id"]=labels.iloc[jj]["Row.names"].to_numpy(); mapped["pathology_label"]=labels.iloc[jj]["label"].to_numpy(); mapped["registration_residual"]=residual
        mapped.to_csv(sec/"count_to_pathology_mapping.tsv.gz",sep="\t",index=False,compression="gzip")
        labels.to_csv(sec/"pathology_labels.tsv.gz",sep="\t",index=False,compression="gzip")
        pd.DataFrame({"spot_id":mapped.count_spot_id,"x":mapped.new_x,"y":mapped.new_y,"pathology_label":mapped.pathology_label}).to_csv(sec/"spots.tsv.gz",sep="\t",index=False,compression="gzip")
        pd.DataFrame({"metric":["n_count_spots","n_selection_spots","n_pathology_spots","mapping_orientation","mapping_median_residual","mapping_p95_residual","mapping_max_residual"],"value":[len(counts),len(selection),len(labels),orientation,float(np.median(residual)),float(np.quantile(residual,.95)),float(np.max(residual))]}).to_csv(sec/"pathology_mapping_summary.tsv",sep="\t",index=False)
        (sec/"section.json").write_text(json.dumps({"cohort":"her2","section":sid,"count_orientation":"spots_by_genes","pathology_mapping":"affine_one_to_one","registration_orientation":orientation},indent=2))
        for src, dst, kind in [(count,sec/"counts.tsv.gz","count_matrix"),(label,sec/"pathology_coordinates.tsv","pathology_coordinates"),(sel_path,sec/"selection.tsv.gz","selection")]:
            manifest.append({"cohort":"her2","section":sid,"kind":kind,"source":str(src),"staged":str(dst),"sha256":sha256(dst)})

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--data-root",type=Path,required=True); ap.add_argument("--output-root",type=Path,required=True); ap.add_argument("--expected-sections",type=int,default=14); a=ap.parse_args()
    out=a.output_root/"inputs"; out.mkdir(parents=True,exist_ok=True); manifest=[]
    stage_atlas(a.data_root,out,manifest); stage_her2(a.data_root,out,manifest)
    sections=sorted((out/"canonical").glob("*/section.json"))
    if len(sections)!=a.expected_sections: raise RuntimeError(f"Expected {a.expected_sections}, found {len(sections)}")
    pd.DataFrame(manifest).to_csv(out/"input_file_manifest.tsv",sep="\t",index=False)
    (out/"PREPARED_OK").write_text("Canonical inputs retained; HER2 pathology mapping intentionally deferred to the explicit mapper.\n")
if __name__=="__main__": main()
