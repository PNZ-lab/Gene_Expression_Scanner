# BeatAMLExplorer

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Literal
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

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

    @staticmethod
    def _sample_cols_from_expr(df: pd.DataFrame) -> list[str]:
        return [c for c in df.columns if c.startswith("BA")]

    def get_clinical_subset(
        self,
        clinical_col: str,
        values: str | list[str],
        match_mode: Literal["exact", "contains"] = "contains",
        sample_col: str = "dbgap_rnaseq_sample",
    ) -> pd.DataFrame:
        """
        Return the clinical rows matching a user-defined subset.

        This is intentionally simple and IDE-friendly: edit the column name,
        values, and match mode in the analysis section below.
        """
        if self.clinical is None:
            self.load_optional_annotations()
        if self.clinical is None:
            raise ValueError("Clinical annotations are not available.")

        if clinical_col not in self.clinical.columns:
            raise ValueError(
                f"'{clinical_col}' not found in clinical table. "
                f"Available columns: {self.clinical.columns.tolist()}"
            )
        if sample_col not in self.clinical.columns:
            raise ValueError(
                f"'{sample_col}' not found in clinical table. "
                f"Available columns: {self.clinical.columns.tolist()}"
            )

        wanted = [values] if isinstance(values, str) else list(values)
        subset = self.clinical.copy()
        series = subset[clinical_col].astype(str)

        if match_mode == "exact":
            mask = series.isin(wanted)
        elif match_mode == "contains":
            upper_vals = [v.upper() for v in wanted]
            mask = series.str.upper().apply(lambda x: any(v in x for v in upper_vals))
        else:
            raise ValueError("match_mode must be 'exact' or 'contains'.")

        return subset.loc[mask].copy()

    def get_clinical_subset_samples(
        self,
        clinical_col: str,
        values: str | list[str],
        match_mode: Literal["exact", "contains"] = "contains",
        sample_col: str = "dbgap_rnaseq_sample",
    ) -> pd.Index:
        subset = self.get_clinical_subset(
            clinical_col=clinical_col,
            values=values,
            match_mode=match_mode,
            sample_col=sample_col,
        )
        return pd.Index(subset[sample_col].dropna().astype(str).unique())

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

    def get_signature_expression(
        self,
        genes_up: list[str],
        genes_down: Optional[list[str]] = None,
        gene_col: str = "display_label",
        zscore: bool = True,
        agg: Literal["first", "mean"] = "mean",
    ) -> pd.Series:
        """
        Compute a signed per-sample signature score.

        Score = mean(z(up genes)) - mean(z(down genes))
        If genes_down is omitted, returns mean(z(up genes)).
        """
        if self.expr is None:
            self.load_expression()

        up_vectors = []
        missing_genes = []
        for gene in genes_up:
            try:
                up_vectors.append(self.get_gene_expression(gene, gene_col=gene_col, agg=agg))
            except ValueError:
                missing_genes.append(gene)

        if not up_vectors:
            raise ValueError(
                "None of the requested up-genes were found in the expression table. "
                f"Missing genes: {missing_genes}"
            )

        up_df = pd.concat(up_vectors, axis=1)
        if zscore:
            up_df = up_df.apply(lambda col: stats.zscore(col, nan_policy="omit"), axis=0, result_type="broadcast")
        up_score = up_df.mean(axis=1)

        if genes_down:
            down_vectors = []
            for gene in genes_down:
                try:
                    down_vectors.append(self.get_gene_expression(gene, gene_col=gene_col, agg=agg))
                except ValueError:
                    missing_genes.append(gene)
            if down_vectors:
                down_df = pd.concat(down_vectors, axis=1)
                if zscore:
                    down_df = down_df.apply(lambda col: stats.zscore(col, nan_policy="omit"), axis=0, result_type="broadcast")
                score = up_score - down_df.mean(axis=1)
            else:
                score = up_score
        else:
            score = up_score

        score = pd.to_numeric(score, errors="coerce")
        score.name = "signature_score"
        score.index.name = "dbgap_rnaseq_sample"
        return score

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

    def _build_feature_drug_frame(
        self,
        feature_vec: pd.Series,
        clinical_col: Optional[str] = None,
        clinical_values: Optional[str | list[str]] = None,
        clinical_match_mode: Literal["exact", "contains"] = "contains",
        sample_ids: Optional[list[str] | pd.Index | set[str]] = None,
    ) -> pd.DataFrame:
        if self.drug_long is None:
            self.build_drug_long()

        df = self.drug_long.copy()
        df[feature_vec.name] = df["dbgap_rnaseq_sample"].map(feature_vec)
        df = df.dropna(subset=[feature_vec.name, "auc"]).copy()

        allowed_samples = None
        if clinical_col is not None and clinical_values is not None:
            allowed_samples = set(
                self.get_clinical_subset_samples(
                    clinical_col=clinical_col,
                    values=clinical_values,
                    match_mode=clinical_match_mode,
                )
            )
        if sample_ids is not None:
            sample_set = set(map(str, sample_ids))
            allowed_samples = sample_set if allowed_samples is None else allowed_samples & sample_set

        if allowed_samples is not None:
            df = df[df["dbgap_rnaseq_sample"].astype(str).isin(allowed_samples)].copy()

        df["sens_score"] = -df["auc"]
        return df

    @staticmethod
    def _corr_with_p(x: pd.Series, y: pd.Series, method: Literal["pearson", "spearman"]) -> tuple[float, float]:
        pair = pd.DataFrame({"x": x, "y": y}).dropna()
        if len(pair) < 3:
            return np.nan, np.nan
        if method == "pearson":
            r, p = stats.pearsonr(pair["x"], pair["y"])
        else:
            r, p = stats.spearmanr(pair["x"], pair["y"])
        return float(r), float(p)

    def feature_drug_sensitivity_scan(
        self,
        feature_vec: pd.Series,
        feature_name: Optional[str] = None,
        min_n: int = 25,
        method: Literal["pearson", "spearman", "both"] = "both",
        split_high_low: bool = True,
        clinical_col: Optional[str] = None,
        clinical_values: Optional[str | list[str]] = None,
        clinical_match_mode: Literal["exact", "contains"] = "contains",
        sample_ids: Optional[list[str] | pd.Index | set[str]] = None,
    ) -> pd.DataFrame:
        """
        Generic feature-to-drug scan used by both gene and signature analyses.

        `feature_vec` should be indexed by `dbgap_rnaseq_sample`.
        """
        feature = feature_vec.copy()
        feature.name = feature_name or feature.name or "feature_value"

        df = self._build_feature_drug_frame(
            feature,
            clinical_col=clinical_col,
            clinical_values=clinical_values,
            clinical_match_mode=clinical_match_mode,
            sample_ids=sample_ids,
        )

        results = []
        for inh, g in df.groupby("inhibitor", dropna=True):
            g2 = g.dropna(subset=[feature.name, "auc"]).copy()
            n = int(g2.shape[0])
            if n < min_n:
                continue

            row = {"feature": feature.name, "inhibitor": inh, "n": n}

            if method in ("pearson", "both"):
                r_auc, p_auc = self._corr_with_p(g2[feature.name], g2["auc"], "pearson")
                r_sens, p_sens = self._corr_with_p(g2[feature.name], g2["sens_score"], "pearson")
                row["r_pearson_expr_auc"] = r_auc
                row["p_pearson_expr_auc"] = p_auc
                row["r_pearson_expr_sens"] = r_sens
                row["p_pearson_expr_sens"] = p_sens

            if method in ("spearman", "both"):
                r_auc_s, p_auc_s = self._corr_with_p(g2[feature.name], g2["auc"], "spearman")
                r_sens_s, p_sens_s = self._corr_with_p(g2[feature.name], g2["sens_score"], "spearman")
                row["r_spearman_expr_auc"] = r_auc_s
                row["p_spearman_expr_auc"] = p_auc_s
                row["r_spearman_expr_sens"] = r_sens_s
                row["p_spearman_expr_sens"] = p_sens_s

            if split_high_low:
                med = float(g2[feature.name].median())
                hi = g2[g2[feature.name] >= med]["auc"]
                lo = g2[g2[feature.name] < med]["auc"]
                row["delta_mean_auc_hi_minus_lo"] = float(hi.mean() - lo.mean())
                if len(hi) >= 2 and len(lo) >= 2:
                    row["p_mannwhitney_hi_vs_lo_auc"] = float(
                        stats.mannwhitneyu(hi, lo, alternative="two-sided").pvalue
                    )
                else:
                    row["p_mannwhitney_hi_vs_lo_auc"] = np.nan

            results.append(row)

        res = pd.DataFrame(results)
        if res.empty:
            return res

        for p_col in [
            "p_spearman_expr_sens",
            "p_pearson_expr_sens",
            "p_mannwhitney_hi_vs_lo_auc",
        ]:
            if p_col in res.columns:
                res[f"fdr_{p_col}"] = self._bh_fdr(res[p_col])

        if "r_spearman_expr_sens" in res.columns:
            res["strength"] = res["r_spearman_expr_sens"].astype(float)
        elif "r_pearson_expr_sens" in res.columns:
            res["strength"] = res["r_pearson_expr_sens"].astype(float)
        elif "r_spearman_expr_auc" in res.columns:
            res["strength"] = -res["r_spearman_expr_auc"].astype(float)
        else:
            res["strength"] = -res["r_pearson_expr_auc"].astype(float)

        return res.sort_values(["strength", "n"], ascending=[False, False]).reset_index(drop=True)

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
        expr_vec = self.get_gene_expression(gene=gene, gene_col=gene_col, agg="mean")
        return self.feature_drug_sensitivity_scan(
            feature_vec=expr_vec,
            feature_name=gene,
            min_n=min_n,
            method=method,
            split_high_low=split_high_low,
        )

    def signature_drug_sensitivity_scan(
        self,
        genes_up: list[str],
        genes_down: Optional[list[str]] = None,
        gene_col: str = "display_label",
        min_n: int = 25,
        method: Literal["pearson", "spearman", "both"] = "both",
        split_high_low: bool = True,
        clinical_col: Optional[str] = None,
        clinical_values: Optional[str | list[str]] = None,
        clinical_match_mode: Literal["exact", "contains"] = "contains",
        sample_ids: Optional[list[str] | pd.Index | set[str]] = None,
        signature_name: str = "signature_score",
    ) -> pd.DataFrame:
        score = self.get_signature_expression(
            genes_up=genes_up,
            genes_down=genes_down,
            gene_col=gene_col,
        )
        score.name = signature_name
        return self.feature_drug_sensitivity_scan(
            feature_vec=score,
            feature_name=signature_name,
            min_n=min_n,
            method=method,
            split_high_low=split_high_low,
            clinical_col=clinical_col,
            clinical_values=clinical_values,
            clinical_match_mode=clinical_match_mode,
            sample_ids=sample_ids,
        )

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
        with optional linear regression and Pearson R + p-value.
        """
        expr_vec = self.get_gene_expression(gene=gene, gene_col=gene_col, agg="mean")
        self.plot_feature_vs_drug(expr_vec, inhibitor, feature_name=gene, label_n=label_n, show_regression=show_regression)

    def plot_feature_vs_drug(
        self,
        feature_vec: pd.Series,
        inhibitor: str,
        feature_name: Optional[str] = None,
        label_n: bool = True,
        show_regression: bool = True,
        clinical_col: Optional[str] = None,
        clinical_values: Optional[str | list[str]] = None,
        clinical_match_mode: Literal["exact", "contains"] = "contains",
        sample_ids: Optional[list[str] | pd.Index | set[str]] = None,
    ) -> None:
        feature = feature_vec.copy()
        feature.name = feature_name or feature.name or "feature_value"
        df = self._build_feature_drug_frame(
            feature,
            clinical_col=clinical_col,
            clinical_values=clinical_values,
            clinical_match_mode=clinical_match_mode,
            sample_ids=sample_ids,
        )
        df = df[df["inhibitor"] == inhibitor].dropna(subset=[feature.name, "auc"]).copy()

        x = df[feature.name].values
        y = df["auc"].values

        plt.figure(figsize=(6, 5), dpi=200)
        plt.scatter(x, y, alpha=0.5)

        title = f"{feature.name} vs {inhibitor}"

        if show_regression and len(df) >= 3:
            slope, intercept = np.polyfit(x, y, 1)
            xfit = np.linspace(x.min(), x.max(), 100)
            yfit = slope * xfit + intercept
            plt.plot(xfit, yfit, color='black')
            r, p = stats.pearsonr(x, y)
            title += f" | R = {r:.2f}, p = {p:.2e}"

        if label_n:
            title += f" (n={len(df)})"

        plt.xlabel(feature.name)
        plt.ylabel("AUC (lower = more sensitive)")
        plt.title(title)
        plt.tight_layout()
        plt.show()

    def plot_signature_vs_drug(
        self,
        genes_up: list[str],
        genes_down: Optional[list[str]],
        inhibitor: str,
        gene_col: str = "display_label",
        signature_name: str = "signature_score",
        clinical_col: Optional[str] = None,
        clinical_values: Optional[str | list[str]] = None,
        clinical_match_mode: Literal["exact", "contains"] = "contains",
        sample_ids: Optional[list[str] | pd.Index | set[str]] = None,
    ) -> None:
        score = self.get_signature_expression(genes_up, genes_down=genes_down, gene_col=gene_col)
        score.name = signature_name
        self.plot_feature_vs_drug(
            score,
            inhibitor=inhibitor,
            feature_name=signature_name,
            clinical_col=clinical_col,
            clinical_values=clinical_values,
            clinical_match_mode=clinical_match_mode,
            sample_ids=sample_ids,
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

    @staticmethod
    def _infer_col(df: pd.DataFrame, candidates: list[str]) -> str:
        cols = set(df.columns)
        for c in candidates:
            if c in cols:
                return c
        raise ValueError(
            f"Could not infer column. Tried: {candidates}. Available: {df.columns.tolist()}"
        )
    
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


def summarize_feature_hi_lo(
    bx,
    feature_vec: pd.Series,
    inhibitor: str,
):
    feature = feature_vec.copy()
    if feature.name is None:
        feature.name = "feature_value"

    df = bx.drug_long[bx.drug_long["inhibitor"] == inhibitor].copy()
    df[feature.name] = df["dbgap_rnaseq_sample"].map(feature)
    df = df.dropna(subset=[feature.name, "auc"])

    med = df[feature.name].median()
    df["group"] = np.where(df[feature.name] >= med, "high", "low")

    return (
        df.groupby("group")["auc"]
        .agg(["count", "mean", "median", "std"])
        .assign(group=lambda x: x.index)
        .reset_index(drop=True)
    )

#%% =============================================================================
# All for one gene
# =============================================================================


gene = 'TYROBP'

bx = BeatAMLExplorer(base_dir=in_path)
bx.load_expression()
bx.load_probit()
bx.load_optional_annotations()

# Scan all drugs for a gene
res = bx.gene_drug_sensitivity_scan(gene, min_n=25, method="both", split_high_low=True)

# Plot a specific hit
top_inhibitor = res.loc[0, "inhibitor"]
bx.plot_gene_vs_drug(gene, top_inhibitor)

# bx.plot_gene_vs_drug(gene, 'GSK-2879552')



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


#%% =============================================================================
# Signature scan for Menin overlap genes
# =============================================================================

# Edit these lists directly in your IDE.
signature_genes_up = [
    "CTSG",
    "CCND1",
]

signature_genes_down = [
    "MEN1",
    "RUNX2",
    "SIPA1L2",
    "PMEPA1",
    "PREX2",
    "NEXN",
    "PI16",
    "RCBTB2",
    "MNS1",
    "HES1",
    "SCN3A",
    "IKZF3",
    "SLA",
    "GUCY1A1",
    "FAXC",
    "TRABD2A",
    "SMAD9",
    "SKIDA1",
]

# Optional clinical subset for KMT2A-rearranged AML.
# Set to None to scan all AML samples.
signature_subset_col = "consensusAMLFusions"
signature_subset_values = ["KMT2A"]
signature_subset_match_mode = "contains"

bx = BeatAMLExplorer(base_dir=in_path)
bx.load_expression()
bx.load_probit()
bx.load_optional_annotations()
bx.build_drug_long()

sig_res = bx.signature_drug_sensitivity_scan(
    genes_up=signature_genes_up,
    genes_down=signature_genes_down,
    min_n=15,
    method="both",
    split_high_low=True,
    clinical_col=signature_subset_col,
    clinical_values=signature_subset_values,
    clinical_match_mode=signature_subset_match_mode,
    signature_name="Menin_overlap_signature",
)

print(sig_res.head(30))

if not sig_res.empty:
    sig_res.to_csv(bx._path("menin_overlap_signature_drug_scan.csv"), index=False)
    top_signature_hit = sig_res.iloc[0]["inhibitor"]
    bx.plot_signature_vs_drug(
        genes_up=signature_genes_up,
        genes_down=signature_genes_down,
        inhibitor=top_signature_hit,
        signature_name="Menin_overlap_signature",
        clinical_col=signature_subset_col,
        clinical_values=signature_subset_values,
        clinical_match_mode=signature_subset_match_mode,
    )
    print("Top signature hit:", top_signature_hit)

