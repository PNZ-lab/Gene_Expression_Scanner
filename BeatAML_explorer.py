# BeatAMLExplorer

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Literal
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

in_path = Path('/Users/kachrist/Desktop/BeatAML')


@dataclass
class BeatAMLExplorer:
    """
    Minimal BeatAML explorer.

    Expects these files in `base_dir`:
      - beataml_waves1to4_norm_exp_dbgap.txt
      - beataml_probit_curve_fits_v4_dbgap.txt
    Optional:
      - beataml_wv1to4_raw_inhibitor_v4_dbgap.txt
      - beataml_drug_families.xlsx
      - beataml_wv1to4_clinical.xlsx
    """

    base_dir: str | Path

    expr_file: str = "beataml_waves1to4_norm_exp_dbgap.txt"
    probit_file: str = "beataml_probit_curve_fits_v4_dbgap.txt"
    inhibitor_file: str = "beataml_wv1to4_raw_inhibitor_v4_dbgap.txt"
    drug_families_file: str = "beataml_drug_families.xlsx"
    clinical_file: str = "beataml_wv1to4_clinical.xlsx"

    expr: Optional[pd.DataFrame] = None          # gene annotations + BAxxxxR columns
    probit: Optional[pd.DataFrame] = None        # drug fits long-ish
    drug_long: Optional[pd.DataFrame] = None     # clean (sample, drug, auc, ic50, ...)
    drug_families: Optional[pd.DataFrame] = None
    clinical: Optional[pd.DataFrame] = None

    def _path(self, name: str) -> Path:
        return Path(self.base_dir) / name

    def load_expression(self) -> pd.DataFrame:
        """
        Loads the normalized expression table.
        Assumes gene identifiers live in columns like stable_id / display_label
        and samples are columns starting with 'BA'.
        """
        p = self._path(self.expr_file)
        df = pd.read_csv(p, sep="\t", low_memory=False)

        # Make it easier to find gene symbols
        # Many BeatAML tables use 'display_label' as gene symbol
        # Keep as-is but normalize column names.
        df.columns = [c.strip() for c in df.columns]

        sample_cols = [c for c in df.columns if c.startswith("BA")]
        if len(sample_cols) == 0:
            raise ValueError(
                "No sample columns found starting with 'BA' in expression file. "
                "Check delimiter and file format."
            )

        self.expr = df
        return df

    def load_probit(self) -> pd.DataFrame:
        """
        Loads the probit curve fits table.
        """
        p = self._path(self.probit_file)
        df = pd.read_csv(p, sep="\t", low_memory=False)
        df.columns = [c.strip() for c in df.columns]
        self.probit = df
        return df

    def load_optional_annotations(self) -> None:
        """
        Loads optional metadata if present.
        """
        # Drug families
        p = self._path(self.drug_families_file)
        if p.exists():
            fam = pd.read_excel(p)
            fam.columns = [c.strip() for c in fam.columns]
            self.drug_families = fam

        # Clinical
        p = self._path(self.clinical_file)
        if p.exists():
            clin = pd.read_excel(p)
            clin.columns = [c.strip() for c in clin.columns]
            self.clinical = clin

    def get_gene_expression(
        self,
        gene: str,
        gene_col: str = "display_label",
        agg: Literal["first", "mean"] = "mean",
    ) -> pd.Series:
        """
        Returns a Series indexed by sample (BAxxxxR), values = expression.

        If multiple rows match the gene (isoforms etc), use mean or first.
        """
        if self.expr is None:
            self.load_expression()

        df = self.expr
        if gene_col not in df.columns:
            raise ValueError(
                f"'{gene_col}' not found in expression table. "
                f"Available columns include: {df.columns[:20].tolist()} ..."
            )

        hits = df[df[gene_col].astype(str).str.upper() == gene.upper()]
        if hits.empty:
            # Try a contains fallback to help debug
            approx = df[df[gene_col].astype(str).str.upper().str.contains(gene.upper(), na=False)].head(10)
            raise ValueError(
                f"Gene '{gene}' not found using gene_col='{gene_col}'. "
                f"Closest matches (first 10): {approx[gene_col].tolist()}"
            )

        sample_cols = [c for c in df.columns if c.startswith("BA")]
        mat = hits[sample_cols]

        if agg == "first":
            vec = mat.iloc[0]
        else:
            vec = mat.mean(axis=0)

        vec.name = gene
        vec.index.name = "dbgap_rnaseq_sample"
        return vec

    def build_drug_long(
        self,
        type_filter: Optional[str] = "single-agent",
        status_filter: Optional[str] = None,
        auc_col_candidates: tuple[str, ...] = ("auc", "AUC"),
        sample_col_candidates: tuple[str, ...] = ("dbgap_rnaseq_sample", "sample", "rnaseq_sample"),
        drug_col_candidates: tuple[str, ...] = ("inhibitor", "drug", "compound"),
    ) -> pd.DataFrame:
        """
        Standardizes probit fits into a clean long table with columns:
          - dbgap_rnaseq_sample
          - inhibitor
          - auc (if present)
          - ic50 (if present)
          - plus any other columns retained
        """
        if self.probit is None:
            self.load_probit()

        df = self.probit.copy()

        # Optional filters if the columns exist
        if type_filter is not None and "type" in df.columns:
            df = df[df["type"] == type_filter]
        if status_filter is not None and "status" in df.columns:
            df = df[df["status"] == status_filter]

        # Identify key columns robustly
        def pick_col(cands: tuple[str, ...]) -> str:
            for c in cands:
                if c in df.columns:
                    return c
            raise ValueError(f"Could not find any of columns: {cands}. Available: {df.columns.tolist()}")

        sample_col = pick_col(sample_col_candidates)
        drug_col = pick_col(drug_col_candidates)

        # AUC column
        auc_col = None
        for c in auc_col_candidates:
            if c in df.columns:
                auc_col = c
                break

        if auc_col is None:
            raise ValueError(f"Could not find AUC column among {auc_col_candidates}. Available: {df.columns.tolist()}")

        out = df.rename(columns={sample_col: "dbgap_rnaseq_sample", drug_col: "inhibitor", auc_col: "auc"}).copy()

        # Drop rows without essentials
        out = out.dropna(subset=["dbgap_rnaseq_sample", "inhibitor", "auc"])

        # Ensure numeric
        out["auc"] = pd.to_numeric(out["auc"], errors="coerce")
        out = out.dropna(subset=["auc"])

        self.drug_long = out
        return out

    def gene_drug_sensitivity_scan(
        self,
        gene: str,
        gene_col: str = "display_label",
        min_n: int = 25,
        method: Literal["pearson", "spearman", "both"] = "both",
        split_high_low: bool = True,
    ) -> pd.DataFrame:
        """
        For each inhibitor, test whether higher gene expression associates with sensitivity.

        Sensitivity is operationalized as:
          - lower AUC (so correlation expr vs AUC should be negative for sensitization)
          - we also report correlation vs (-AUC) as "sens_score"

        Returns a dataframe sorted by strongest sensitization evidence.
        """
        if self.drug_long is None:
            self.build_drug_long()

        expr_vec = self.get_gene_expression(gene=gene, gene_col=gene_col, agg="mean")

        df = self.drug_long.copy()
        df["gene_expr"] = df["dbgap_rnaseq_sample"].map(expr_vec)

        df = df.dropna(subset=["gene_expr"])
        df["sens_score"] = -df["auc"]

        results = []

        for inh, g in df.groupby("inhibitor", dropna=True):
            g2 = g.dropna(subset=["gene_expr", "auc"])
            n = int(g2.shape[0])
            if n < min_n:
                continue

            row = {"inhibitor": inh, "n": n}

            # Correlations
            if method in ("pearson", "both"):
                r_auc = np.corrcoef(g2["gene_expr"], g2["auc"])[0, 1]
                r_sens = np.corrcoef(g2["gene_expr"], g2["sens_score"])[0, 1]
                row["r_pearson_expr_auc"] = r_auc
                row["r_pearson_expr_sens"] = r_sens

            if method in ("spearman", "both"):
                # Spearman using rank correlation via pandas
                r_auc_s = g2["gene_expr"].rank().corr(g2["auc"].rank())
                r_sens_s = g2["gene_expr"].rank().corr(g2["sens_score"].rank())
                row["r_spearman_expr_auc"] = r_auc_s
                row["r_spearman_expr_sens"] = r_sens_s

            # High vs low split effect size
            if split_high_low:
                med = float(g2["gene_expr"].median())
                hi = g2[g2["gene_expr"] >= med]["auc"]
                lo = g2[g2["gene_expr"] < med]["auc"]

                # Simple effect size: delta mean AUC (hi - lo)
                # Negative delta means: high expr has lower AUC (sensitization)
                row["delta_mean_auc_hi_minus_lo"] = float(hi.mean() - lo.mean())

            results.append(row)

        res = pd.DataFrame(results)

        if res.empty:
            return res

        # Prefer sorting by sensitization score:
        # r_expr_sens high is good (since sens_score=-AUC)
        sort_cols = [c for c in ["r_spearman_expr_sens", "r_pearson_expr_sens"] if c in res.columns]
        if sort_cols:
            res = res.sort_values(sort_cols[0], ascending=False)
        else:
            # fallback: negative correlation with AUC indicates sensitization
            sort_cols = [c for c in ["r_spearman_expr_auc", "r_pearson_expr_auc"] if c in res.columns]
            res = res.sort_values(sort_cols[0], ascending=True)

        return res.reset_index(drop=True)

    def plot_gene_vs_drug(
        self,
        gene: str,
        inhibitor: str,
        gene_col: str = "display_label",
        label_n: bool = True,
        show_regression: bool = True,
    ) -> None:
        """
        Scatter plot: gene expression vs AUC for one drug,
        with optional linear regression and Pearson R.
        """
        if self.drug_long is None:
            self.build_drug_long()
    
        expr_vec = self.get_gene_expression(gene=gene, gene_col=gene_col, agg="mean")
    
        df = self.drug_long[self.drug_long["inhibitor"] == inhibitor].copy()
        df["gene_expr"] = df["dbgap_rnaseq_sample"].map(expr_vec)
        df = df.dropna(subset=["gene_expr", "auc"])
    
        x = df["gene_expr"].values
        y = df["auc"].values
    
        plt.figure(figsize=(6, 5), dpi=200)
        plt.scatter(x, y, alpha=0.5)
    
        title = f"{gene} vs {inhibitor}"
    
        if show_regression and len(df) >= 3:
            # Linear regression
            slope, intercept = np.polyfit(x, y, 1)
            xfit = np.linspace(x.min(), x.max(), 100)
            yfit = slope * xfit + intercept
            plt.plot(xfit, yfit, color='black')
    
            # Pearson correlation
            r = np.corrcoef(x, y)[0, 1]
            title += f" | R = {r:.2f}"
    
        if label_n:
            title += f" (n={len(df)})"
    
        plt.xlabel(f"{gene} expression")
        plt.ylabel("AUC (lower = more sensitive)")
        plt.title(title)
        plt.tight_layout()
        plt.show()
        
    @staticmethod
    def _bh_fdr(pvals: pd.Series) -> pd.Series:
        """Benjamini-Hochberg FDR for a pandas Series of p-values."""
        p = pvals.astype(float).values
        n = np.sum(~np.isnan(p))
        out = np.full_like(p, np.nan, dtype=float)
        if n == 0:
            return pd.Series(out, index=pvals.index)

        idx = np.where(~np.isnan(p))[0]
        pv = p[idx]
        order = np.argsort(pv)
        pv_sorted = pv[order]
        ranks = np.arange(1, len(pv_sorted) + 1)
        q_sorted = pv_sorted * len(pv_sorted) / ranks
        q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
        q = np.empty_like(pv_sorted)
        q[order] = np.clip(q_sorted, 0, 1)
        out[idx] = q
        return pd.Series(out, index=pvals.index)

    @staticmethod
    def _infer_col(df: pd.DataFrame, candidates: list[str]) -> str:
        cols = set(df.columns)
        for c in candidates:
            if c in cols:
                return c
        raise ValueError(
            f"Could not infer column. Tried: {candidates}. Available: {df.columns.tolist()}"
        )
    
    @staticmethod
    def _bh_fdr(pvals: pd.Series) -> pd.Series:
        """Benjamini-Hochberg FDR for a pandas Series of p-values."""
        p = pvals.astype(float).values
        n = np.sum(~np.isnan(p))
        out = np.full_like(p, np.nan, dtype=float)
        if n == 0:
            return pd.Series(out, index=pvals.index)

        idx = np.where(~np.isnan(p))[0]
        pv = p[idx]
        order = np.argsort(pv)
        pv_sorted = pv[order]
        ranks = np.arange(1, len(pv_sorted) + 1)
        q_sorted = pv_sorted * len(pv_sorted) / ranks
        q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
        q = np.empty_like(pv_sorted)
        q[order] = np.clip(q_sorted, 0, 1)
        out[idx] = q
        return pd.Series(out, index=pvals.index)


def summarize_hi_lo(
    bx,
    gene: str,
    inhibitor: str,
    gene_col: str = "display_label",
):
    df = bx.drug_long[bx.drug_long["inhibitor"] == inhibitor].copy()
    expr = bx.get_gene_expression(gene, gene_col=gene_col)

    df["expr"] = df["dbgap_rnaseq_sample"].map(expr)
    df = df.dropna(subset=["expr", "auc"])

    med = df["expr"].median()
    df["group"] = np.where(df["expr"] >= med, "high", "low")

    summary = (
        df.groupby("group")["auc"]
        .agg(["count", "mean", "median", "std"])
        .assign(group=lambda x: x.index)
        .reset_index(drop=True)
    )

    return summary

#%% =============================================================================
# All for one gene
# =============================================================================

gene = 'EZH2'

bx = BeatAMLExplorer(base_dir=in_path)
bx.load_expression()
bx.load_probit()
bx.load_optional_annotations()

# Scan all drugs for a gene
res = bx.gene_drug_sensitivity_scan(gene, min_n=25, method="both", split_high_low=True)

# Plot a specific hit
top_inhibitor = res.loc[0, "inhibitor"]
bx.plot_gene_vs_drug(gene, top_inhibitor)

print(gene)
print(res)
print(top_inhibitor)

summarize_hi_lo(bx, gene, top_inhibitor)


expr = bx.expr if bx.expr is not None else bx.load_expression()



#%% =============================================================================
# Calculate all genes to all drugs
# =============================================================================

# BeatAMLExplorer
import sys
from tqdm import tqdm

# -----------------------
# Setup
# -----------------------
bx = BeatAMLExplorer(base_dir=in_path)
bx.load_expression()
bx.load_probit()
bx.load_optional_annotations()
bx.build_drug_long()  # build once

expr = bx.expr

# Gene list: use display_label (symbols). Drop blanks/NA.
gene_col = "display_label"
gene_list = (
    expr[gene_col]
    .dropna()
    .astype(str)
    .str.strip()
)
gene_list = gene_list[gene_list != ""].unique().tolist()

print("N genes (unique symbols):", len(gene_list))
print("Example:", gene_list[:10])

# -----------------------
# Scan parameters
# -----------------------
min_n = 25
method = "both"
split_high_low = True

# Keep only top K drugs per gene (prevents huge memory use)
TOPK_PER_GENE = 3

all_top_hits = []
errors = []

pbar = tqdm(gene_list, desc="Scanning genes", file=sys.stdout)
for gene in pbar:
    pbar.set_postfix_str(gene)
    try:
        res = bx.gene_drug_sensitivity_scan(
            gene=gene,
            gene_col=gene_col,
            min_n=min_n,
            method=method,
            split_high_low=split_high_low
        )
        if res is None or res.empty:
            continue

        # Add gene label
        res = res.copy()
        res.insert(0, "gene", gene)

        # Choose a single "strength" column to rank within gene
        # Prefer spearman on sens_score, else pearson on sens_score.
        if "r_spearman_expr_sens" in res.columns:
            res["strength"] = res["r_spearman_expr_sens"].astype(float)
        elif "r_pearson_expr_sens" in res.columns:
            res["strength"] = res["r_pearson_expr_sens"].astype(float)
        else:
            # fallback: negative corr with AUC means sensitization
            if "r_spearman_expr_auc" in res.columns:
                res["strength"] = (-res["r_spearman_expr_auc"].astype(float))
            else:
                res["strength"] = (-res["r_pearson_expr_auc"].astype(float))

        # Keep only the strongest TOPK per gene
        res = res.sort_values("strength", ascending=False).head(TOPK_PER_GENE)

        all_top_hits.append(res)

    except Exception as e:
        errors.append({"gene": gene, "error": repr(e)})

hits_df = pd.concat(all_top_hits, ignore_index=True) if all_top_hits else pd.DataFrame()
errors_df = pd.DataFrame(errors)

print("Genes scanned:", len(gene_list))
print("Genes with >=1 hit:", hits_df["gene"].nunique() if not hits_df.empty else 0)
print("Errors:", len(errors_df))

# -----------------------
# Global ranking: "most sensitivity-determining"
# -----------------------
# This ranks across all (gene, inhibitor) pairs.
# You can tighten thresholds to keep only very strong signals.
if not hits_df.empty:
    hits_df = hits_df.sort_values(["strength", "n"], ascending=[False, False]).reset_index(drop=True)

    # Optional: filter for "highly correlated"
    # e.g. keep abs(corr) >= 0.4
    CORR_MIN = 0.4
    hits_df = hits_df[hits_df["strength"].abs() >= CORR_MIN].reset_index(drop=True)

print("\nTop 30 strongest gene–drug sensitivity associations:")
print(hits_df.head(30))

# Save results
hits_df.to_csv(bx._path("gene_drug_scan_top_hits.csv"), index=False)
errors_df.to_csv(bx._path("gene_drug_scan_errors.csv"), index=False)

print("\nSaved:")
print(" - gene_drug_scan_top_hits.csv")
print(" - gene_drug_scan_errors.csv")

#%%
TOP_N = 30
top_hits = hits_df.head(TOP_N)
for i, row in top_hits.iterrows():
    print(f"{i+1}/{TOP_N}: {row['gene']} vs {row['inhibitor']} "
          f"(strength={row['strength']:.2f}, n={row['n']})")
    bx.plot_gene_vs_drug(
        gene=row["gene"],
        inhibitor=row["inhibitor"],
    )

