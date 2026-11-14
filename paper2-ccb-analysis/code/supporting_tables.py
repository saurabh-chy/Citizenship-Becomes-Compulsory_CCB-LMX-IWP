"""
Paper 2 - supporting tables (descriptives, reliability, curvilinear models,
H5/H6 detail, regression diagnostics, VIFs and sensitivity checks).

These analyses come unchanged from the thesis Chapter 4 pipeline (Appendix D.2).
They do not depend on the bootstrap, so they are unaffected by the later
mediator-unit rescaling. H4 and H7 bootstrap estimates and all paper figures
are produced by paper2_analysis.py instead.

Usage:
    python supporting_tables.py --data path/to/CCB_Two_Phase_Data.csv --out outputs/supporting_tables
"""
import argparse, re, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import FactorAnalysis
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan, linear_reset
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--data", required=True)
ap.add_argument("--out", default="outputs/supporting_tables")
args = ap.parse_args()
DATA_PATH = Path(args.data)
TABLE_DIR = Path(args.out); TABLE_DIR.mkdir(parents=True, exist_ok=True)
SEED = 20260825

# Load the dataset
df = pd.read_excel(DATA_PATH, sheet_name="Original_Data") if DATA_PATH.suffix.lower() in (".xlsx", ".xls") else pd.read_csv(DATA_PATH)
if len(df) != 391:
    raise ValueError(f"Expected 391 matched records, found {len(df)}")
# 3. Variable map and reproducible composite scoring
 
SCALES = {"CCB":5,"LMX1":7,"LMX2":7,"TP":5,"CP":8,"CWB":5,
          "SE":10,"OR":10,"DU":10,"AS":10,"SD":10,"CA":10}
CORE_SCALES={k:SCALES[k] for k in ["CCB","LMX1","LMX2","TP","CP","CWB"]}
FACETS=["SE","OR","DU","AS","SD","CA"]
CONTROLS=["Gender_Female","Marital_Status_Married","Education_Level","Total_Experience_Years","Tenure_Org_Months","Tenure_Supervisor_Months"]
MIN_ITEM_PROPORTION=.80
def items(prefix,n): return [f"{prefix}_{i}" for i in range(1,n+1)]
required=[x for p,n in SCALES.items() for x in items(p,n)]+["Age"]+CONTROLS
missing=sorted(set(required)-set(df.columns))
if missing: raise KeyError(f"Required columns missing: {missing}")
score_audit=[]
for prefix,nitem in SCALES.items():
    cols=items(prefix,nitem); valid=df[cols].notna().sum(axis=1); minimum=int(np.ceil(MIN_ITEM_PROPORTION*nitem))
    df[f"{prefix}_mean"]=df[cols].mean(axis=1).where(valid>=minimum)
    score_audit.append({"Construct":prefix,"Items":nitem,"Minimum_valid_items":minimum,"Scores_set_missing":int((valid<minimum).sum())})
pd.DataFrame(score_audit).to_csv(TABLE_DIR/"Table_Scoring_Audit.csv",index=False)
df["CWB_reverse"]=6-df["CWB_mean"]
df["CONSC_mean"]=df[[f"{x}_mean" for x in FACETS]].mean(axis=1)
df["IWP_mean"]=(df["TP_mean"]+df["CP_mean"]+df["CWB_reverse"])/3
levels=sorted(df.Education_Level.dropna().unique())
for level in levels[1:]: df[f"Education_{int(level)}"]=(df.Education_Level==level).astype(int)
def zscore(x): return (x-x.mean())/x.std(ddof=0)
for v in ["CCB_mean","LMX1_mean","LMX2_mean","TP_mean","CP_mean","CWB_mean","CONSC_mean","IWP_mean"]: df[f"z_{v}"]=zscore(df[v])
df["CCB_sq"]=df.z_CCB_mean**2; df["CCB_x_CONSC"]=df.z_CCB_mean*df.z_CONSC_mean; df["CCBsq_x_CONSC"]=df.CCB_sq*df.z_CONSC_mean
for facet in FACETS:
    df[f"z_{facet}_mean"]=zscore(df[f"{facet}_mean"]); df[f"CCBsq_x_{facet}"]=df.CCB_sq*df[f"z_{facet}_mean"]
 
# 4. Data-quality, eligibility, and distribution checks
 
all_items = [c for prefix, n_items in SCALES.items() for c in items(prefix, n_items)]
checks = pd.DataFrame({
    "Check": ["Rows", "Columns", "Duplicate rows", "Missing cells", "Item values below 1", "Item values above 5",
              "Supervisor tenure > organization tenure", "Organization tenure > total experience (months)",
              "Age < total experience + 15"],
    "Result": [len(df), len(df.columns), int(df.duplicated().sum()), int(df.isna().sum().sum()),
               int((df[all_items] < 1).sum().sum()), int((df[all_items] > 5).sum().sum()),
               int((df["Tenure_Supervisor_Months"] > df["Tenure_Org_Months"]).sum()),
               int((df["Tenure_Org_Months"] > df["Total_Experience_Years"] * 12).sum()),
               int((df["Age"] < df["Total_Experience_Years"] + 15).sum())]
})
checks.to_csv(TABLE_DIR / "Table_S1_data_quality.csv", index=False)
 
# 5. Descriptive statistics and correlations: Table 1 / Table S2 / Table S5a
 
main_vars = ["CCB_mean", "TP_mean", "CP_mean", "CWB_mean", "LMX1_mean", "LMX2_mean", "CONSC_mean", "IWP_mean"]
labels = {"CCB_mean":"CCB", "TP_mean":"Task performance", "CP_mean":"Contextual performance",
          "CWB_mean":"Counterproductive work behavior", "LMX1_mean":"LMX Time 1",
          "LMX2_mean":"LMX Time 2", "CONSC_mean":"Conscientiousness", "IWP_mean":"IWP composite"}
 
desc = df[main_vars].agg(["count", "mean", "std", "min", "max"]).T.reset_index().rename(columns={"index":"Variable", "count":"N", "mean":"M", "std":"SD", "min":"Min", "max":"Max"})
desc["Variable"] = desc["Variable"].map(labels)
corr = df[main_vars].corr()
corr.index = [labels[x] for x in corr.index]; corr.columns = [labels[x] for x in corr.columns]
table1 = desc.copy()
for c in corr.columns:
    table1[c] = corr[c].to_numpy()
table1.to_csv(TABLE_DIR / "Table_1_descriptives_correlations.csv", index=False)
 
demographics = df[["Age", "Gender_Female", "Marital_Status_Married", "Education_Level", "Total_Experience_Years", "Tenure_Org_Months", "Tenure_Supervisor_Months"]].agg(["count", "mean", "std", "min", "median", "max"]).T
demographics.to_csv(TABLE_DIR / "Table_S2_sample_descriptives.csv")
 
# 6. Reliability analysis
# Python reports alpha, one-factor omega sensitivity, corrected item-total correlations, and alpha if deleted. CFA-based loadings, CR, AVE, Fornell-Larcker, HTMT, invariance, and SEM come only from the companion R/lavaan pipeline.
 
def cronbach_alpha(x):
    x=np.asarray(x,float); k=x.shape[1]
    return k/(k-1)*(1-x.var(axis=0,ddof=1).sum()/x.sum(axis=1).var(ddof=1))
def omega_one_factor(x):
    z=StandardScaler().fit_transform(np.asarray(x,float)); fa=FactorAnalysis(n_components=1,random_state=SEED).fit(z)
    lam=np.abs(fa.components_[0]); uniq=np.maximum(fa.noise_variance_,1e-8)
    return lam.sum()**2/(lam.sum()**2+uniq.sum())
rel_rows=[]; item_rows=[]
for prefix,nitem in SCALES.items():
    X=df[items(prefix,nitem)].dropna(); total=X.sum(axis=1)
    rel_rows.append({"Construct":prefix,"Items":nitem,"N":len(X),"Cronbach_alpha":cronbach_alpha(X),"Omega_one_factor_sensitivity":omega_one_factor(X),"Mean_interitem_correlation":X.corr().values[np.triu_indices(nitem,1)].mean()})
    for item in X:
        rest=total-X[item]
        item_rows.append({"Construct":prefix,"Item":item,"Corrected_item_total_r":X[item].corr(rest),"Alpha_if_deleted":cronbach_alpha(X.drop(columns=item))})
reliability=pd.DataFrame(rel_rows); item_diagnostics=pd.DataFrame(item_rows)
reliability.to_csv(TABLE_DIR/"Table_2_reliability.csv",index=False); item_diagnostics.to_csv(TABLE_DIR/"Table_S4_item_reliability_diagnostics.csv",index=False)
 
# 7. HC3-robust regression helpers
 
def fit_hc3(outcome, predictors, label=None):
    data = df[[outcome] + predictors].dropna().copy()
    X = sm.add_constant(data[predictors], has_constant="add")
    model = sm.OLS(data[outcome], X).fit(cov_type="HC3")
    table = pd.DataFrame({"Predictor": model.params.index, "b": model.params.values, "HC3_SE": model.bse.values,
                          "z_or_t": model.tvalues.values, "p": model.pvalues.values})
    table["R2"] = model.rsquared; table["Adj_R2"] = model.rsquared_adj; table["N"] = int(model.nobs)
    if label:
        table.insert(0, "Outcome", label)
    return model, table
 
def significance(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""
 
# Controls are standardized only for comparability, while binary coding is retained.
for col in ["Total_Experience_Years", "Tenure_Org_Months", "Tenure_Supervisor_Months"]:
    df[f"z_{col}"] = zscore(df[col])
REG_CONTROLS = ["Gender_Female", "Marital_Status_Married", "Education_2", "Education_3", "z_Total_Experience_Years", "z_Tenure_Org_Months", "z_Tenure_Supervisor_Months"]
# Models use already-standardized focal variables.
 
# 8. Curvilinear performance models: Table 4 and Table 5
 
outcomes = {"TP_mean":"Task performance", "CP_mean":"Contextual performance", "CWB_mean":"Counterproductive work behavior", "IWP_mean":"Operationalized IWP"}
quadratic_models, regression_tables, curve_summary = {}, [], []
for outcome, label in outcomes.items():
    predictors = REG_CONTROLS + ["z_CCB_mean", "CCB_sq"]
    model, tab = fit_hc3(outcome, predictors, label)
    quadratic_models[outcome] = model
    regression_tables.append(tab)
    b1, b2 = model.params["z_CCB_mean"], model.params["CCB_sq"]
    turning = -b1 / (2 * b2) if abs(b2) > 1e-12 else np.nan
    xlo, xhi = df["z_CCB_mean"].min(), df["z_CCB_mean"].max()
    slope_lo, slope_hi = b1 + 2*b2*xlo, b1 + 2*b2*xhi
    pattern = "Monotonic decline" if slope_lo < 0 and slope_hi < 0 else ("Inverted-U" if b2 < 0 else "U-shape")
    curve_summary.append({"Outcome":label, "Linear_b":b1, "Linear_p":model.pvalues["z_CCB_mean"],
                          "Quadratic_b":b2, "Quadratic_p":model.pvalues["CCB_sq"],
                          "Turning_point_z":turning, "Slope_low_CCB":slope_lo, "Slope_high_CCB":slope_hi,
                          "Pattern":pattern, "R2":model.rsquared})
 
table4 = pd.concat(regression_tables, ignore_index=True)
table5 = pd.DataFrame(curve_summary)
table4.to_csv(TABLE_DIR / "Table_4_HC3_quadratic_models.csv", index=False)
table5.to_csv(TABLE_DIR / "Table_5_curvilinear_summary.csv", index=False)
 
# 11. Final confirmatory analyses for H1 to H7
# These models are the authoritative Python hypothesis tests. HC3 inference is used throughout. Ordinal CFA, measurement invariance, CFA-based validity, and latent SEM remain in the companion R/lavaan pipeline.
 
# Helpers for confidence intervals, model summaries, nested R-squared changes, and VIF.
def coef_row(model, hypothesis, outcome, model_name, term):
    ci = model.conf_int().loc[term]
    return {"Hypothesis":hypothesis,"Outcome":outcome,"Model":model_name,"Predictor":term,
            "b":model.params[term],"HC3_SE":model.bse[term],"Robust_z":model.tvalues[term],
            "p":model.pvalues[term],"CI_low":ci.iloc[0],"CI_high":ci.iloc[1],
            "N":int(model.nobs),"R2":model.rsquared,"Adjusted_R2":model.rsquared_adj}
 
def nested_change(reduced, full):
    n, q = int(full.nobs), int(full.df_model - reduced.df_model)
    delta = full.rsquared - reduced.rsquared
    f = (delta/q) / ((1-full.rsquared)/(n-full.df_model-1)) if q > 0 else np.nan
    return delta, f, stats.f.sf(f, q, n-full.df_model-1) if q > 0 else np.nan
 
def model_vif(model):
    X = pd.DataFrame(model.model.exog, columns=model.model.exog_names)
    rows=[]
    for col in [c for c in X if c != "const"]:
        others=[c for c in X if c not in ("const",col)]
        r2=sm.OLS(X[col],sm.add_constant(X[others],has_constant="add")).fit().rsquared
        rows.append({"Model_outcome":model.model.endog_names,"Predictor":col,"VIF":1/(1-r2)})
    return pd.DataFrame(rows)
 
 
 
def bca_interval(draws,theta,jack,alpha=.05):
    draws=np.asarray(draws,float); draws=draws[np.isfinite(draws)]; jack=np.asarray(jack,float); jack=jack[np.isfinite(jack)]
    eps=1/(2*len(draws)); prop=np.clip(np.mean(draws<theta),eps,1-eps); z0=stats.norm.ppf(prop)
    jm=jack.mean(); dif=jm-jack; den=6*(np.sum(dif**2)**1.5); acc=np.sum(dif**3)/den if den>0 else 0
    za=stats.norm.ppf([alpha/2,1-alpha/2]); adj=stats.norm.cdf(z0+(z0+za)/(1-acc*(z0+za)))
    return np.quantile(draws,np.clip(adj,0,1)),acc,z0
 
# H1 and H5: baseline-adjusted later LMX
 
# H1 is tested separately from H5.
h1_terms = REG_CONTROLS + ["z_LMX1_mean", "z_CCB_mean"]
h1_model, _ = fit_hc3("LMX2_mean", h1_terms, "H1")
h5_reduced_terms = REG_CONTROLS + ["z_LMX1_mean", "z_CCB_mean", "z_CONSC_mean"]
h5_full_terms = h5_reduced_terms + ["CCB_x_CONSC"]
h5_reduced, _ = fit_hc3("LMX2_mean", h5_reduced_terms, "H5 reduced")
h5_model, _ = fit_hc3("LMX2_mean", h5_full_terms, "H5 full")
dR2_h5, F_h5, pchg_h5 = nested_change(h5_reduced,h5_model)
h15 = pd.DataFrame([
    coef_row(h1_model,"H1","LMX Time 2","Baseline-adjusted", "z_CCB_mean"),
    coef_row(h5_model,"H5","LMX Time 2","Moderation", "CCB_x_CONSC")])
h15["Delta_R2"]=[np.nan,dR2_h5]; h15["Change_F"]=[np.nan,F_h5]; h15["Change_p"]=[np.nan,pchg_h5]
h15.to_csv(TABLE_DIR/"Table_H1_H5_LMX_models.csv",index=False)
 
# Conditional H5 slopes and predicted lines.
h5_slopes=[]
for w in [-1.,0.,1.]:
    slope=h5_model.params["z_CCB_mean"]+w*h5_model.params["CCB_x_CONSC"]
    V=h5_model.cov_params(); var=V.loc["z_CCB_mean","z_CCB_mean"]+w*w*V.loc["CCB_x_CONSC","CCB_x_CONSC"]+2*w*V.loc["z_CCB_mean","CCB_x_CONSC"]
    se=np.sqrt(var); z=slope/se
    h5_slopes.append({"Conscientiousness_SD":w,"Slope":slope,"HC3_SE":se,"CI_low":slope-1.96*se,"CI_high":slope+1.96*se,"p":2*stats.norm.sf(abs(z))})
pd.DataFrame(h5_slopes).to_csv(TABLE_DIR/"Table_H5_conditional_slopes.csv",index=False)
# H2, H3, and H4: performance, curvilinearity, and mediation
 
# H2 dedicated association model.
h2_terms=REG_CONTROLS+["z_LMX1_mean","z_CCB_mean","z_LMX2_mean"]
h2_model,_=fit_hc3("IWP_mean",h2_terms,"H2")
h2_table=pd.DataFrame([coef_row(h2_model,"H2","IWP","Adjusted association","z_LMX2_mean")])
h2_table.to_csv(TABLE_DIR/"Table_H2_LMX_IWP.csv",index=False)
 
# H3 hierarchical model comparisons and strict boundary tests.
def boundary_slope(model,x,positive=True):
    b1,b2=model.params["z_CCB_mean"],model.params["CCB_sq"]; V=model.cov_params(); slope=b1+2*x*b2
    se=np.sqrt(V.loc["z_CCB_mean","z_CCB_mean"]+4*x*x*V.loc["CCB_sq","CCB_sq"]+4*x*V.loc["z_CCB_mean","CCB_sq"])
    z=slope/se; return slope,se,z,(stats.norm.sf(z) if positive else stats.norm.cdf(z))
h3_rows=[]
for y,label in outcomes.items():
    linear,_=fit_hc3(y,REG_CONTROLS+["z_CCB_mean"],label)
    quad=quadratic_models[y]; dr2,fc,pc=nested_change(linear,quad)
    b2=quad.params["CCB_sq"]; concave=b2<0; lo,hi=df.z_CCB_mean.min(),df.z_CCB_mean.max()
    slo=boundary_slope(quad,lo,positive=concave); shi=boundary_slope(quad,hi,positive=not concave)
    tp=-quad.params["z_CCB_mean"]/(2*b2); strict=max(slo[3],shi[3])<.05 and lo<tp<hi and quad.pvalues["CCB_sq"]<.05
    h3_rows.append({**coef_row(quad,"H3",label,"Quadratic","CCB_sq"),"Delta_R2":dr2,"Change_F":fc,"Change_p":pc,
        "Slope_lower":slo[0],"p_lower_one_sided":slo[3],"Slope_upper":shi[0],"p_upper_one_sided":shi[3],
        "Turning_point_z":tp,"Turning_inside_range":lo<tp<hi,"Strict_shape_supported":strict,
        "Curve_classification":("Inverted-U" if concave else "U-shape") if strict else "Strict U-shape not supported"})
h3_table=pd.DataFrame(h3_rows); h3_table.to_csv(TABLE_DIR/"Table_H3_hierarchical_strict_curvilinear_tests.csv",index=False)
 
# H6 and H7: quadratic moderation and moderated mediation
 
# H6 hierarchical interaction tests, robust joint Wald test, conditional effects, and turning points.
h6_reduced_terms=REG_CONTROLS+["z_CCB_mean","CCB_sq","z_CONSC_mean"]
h6_full_terms=h6_reduced_terms+["CCB_x_CONSC","CCBsq_x_CONSC"]
h6_reduced,_=fit_hc3("IWP_mean",h6_reduced_terms,"H6 reduced"); h6_model,_=fit_hc3("IWP_mean",h6_full_terms,"H6 full")
dr2,fchg,pchg=nested_change(h6_reduced,h6_model)
names_h6=h6_model.model.exog_names; R=np.zeros((2,len(names_h6))); R[0,names_h6.index("CCB_x_CONSC")]=1; R[1,names_h6.index("CCBsq_x_CONSC")]=1
h6_wald=h6_model.wald_test(R,use_f=False,scalar=True); h6_wald_chi2=float(h6_wald.statistic); h6_wald_p=float(h6_wald.pvalue)
h6_rows=[]
for term in ["CCB_x_CONSC","CCBsq_x_CONSC"]:
    row=coef_row(h6_model,"H6","IWP","Quadratic moderation",term); row.update({"Delta_R2":dr2,"Conventional_change_F":fchg,"Conventional_change_p":pchg,"HC3_joint_Wald_chi2":h6_wald_chi2,"HC3_joint_Wald_df":2,"HC3_joint_Wald_p":h6_wald_p}); h6_rows.append(row)
 
V=h6_model.cov_params(); params=h6_model.params; lo_x,hi_x=df.z_CCB_mean.min(),df.z_CCB_mean.max(); h6_cond=[]
def linear_combo_h6(weights):
    names=list(weights); estimate=sum(weights[k]*params[k] for k in names); variance=sum(weights[a]*weights[b]*V.loc[a,b] for a in names for b in names); se=np.sqrt(max(variance,0)); z=estimate/se if se>0 else np.nan
    return estimate,se,estimate-1.96*se,estimate+1.96*se,2*stats.norm.sf(abs(z))
for w in [-1.,0.,1.]:
    lin=linear_combo_h6({"z_CCB_mean":1,"CCB_x_CONSC":w}); quad=linear_combo_h6({"CCB_sq":1,"CCBsq_x_CONSC":w})
    L,Q=lin[0],quad[0]; turning=-L/(2*Q) if abs(Q)>1e-12 else np.nan
    if np.isfinite(turning):
        g={"z_CCB_mean":-1/(2*Q),"CCB_x_CONSC":-w/(2*Q),"CCB_sq":L/(2*Q**2),"CCBsq_x_CONSC":w*L/(2*Q**2)}; tse=linear_combo_h6(g)[1]
    else: tse=np.nan
    low_s=linear_combo_h6({"z_CCB_mean":1,"CCB_x_CONSC":w,"CCB_sq":2*lo_x,"CCBsq_x_CONSC":2*lo_x*w})
    high_s=linear_combo_h6({"z_CCB_mean":1,"CCB_x_CONSC":w,"CCB_sq":2*hi_x,"CCBsq_x_CONSC":2*hi_x*w})
    h6_cond.append({"Conscientiousness_SD":w,"Conditional_linear":L,"Linear_SE":lin[1],"Linear_CI_low":lin[2],"Linear_CI_high":lin[3],"Linear_p":lin[4],"Conditional_quadratic":Q,"Quadratic_SE":quad[1],"Quadratic_CI_low":quad[2],"Quadratic_CI_high":quad[3],"Quadratic_p":quad[4],"Turning_point_z":turning,"Turning_SE_delta":tse,"Turning_CI_low":turning-1.96*tse if np.isfinite(tse) else np.nan,"Turning_CI_high":turning+1.96*tse if np.isfinite(tse) else np.nan,"Turning_inside_observed_range":bool(lo_x<turning<hi_x) if np.isfinite(turning) else False,"Slope_at_lower_boundary":low_s[0],"Lower_slope_SE":low_s[1],"Lower_slope_p_two_sided":low_s[4],"Slope_at_upper_boundary":high_s[0],"Upper_slope_SE":high_s[1],"Upper_slope_p_two_sided":high_s[4]})
pd.DataFrame(h6_rows).to_csv(TABLE_DIR/"Table_H6_quadratic_moderation.csv",index=False); h6_cond=pd.DataFrame(h6_cond); h6_cond.to_csv(TABLE_DIR/"Table_H6_conditional_curves.csv",index=False)
 
# H7 point-estimate models (needed for diagnostics; bootstrap is in paper2_analysis.py)
a_terms=REG_CONTROLS+["z_LMX1_mean","z_CCB_mean","z_CONSC_mean","CCB_x_CONSC"]
b_terms=REG_CONTROLS+["z_LMX1_mean","z_CCB_mean","z_CONSC_mean","z_LMX2_mean"]
a7,_=fit_hc3("LMX2_mean",a_terms,"H7 first stage"); b7,_=fit_hc3("IWP_mean",b_terms,"H7 outcome")
# Diagnostics, sensitivity analyses, final workbook, and audit
 
# Diagnostics and model-specific VIFs for the principal models.
principal={"H1":h1_model,"H2":h2_model,"H3_IWP":quadratic_models["IWP_mean"],"H5":h5_model,"H6":h6_model,"H7_a":a7,"H7_b":b7}
diag_rows=[]; vif_rows=[]
for name,m in principal.items():
    infl=m.get_influence(); cook=infl.cooks_distance[0]; lev=infl.hat_matrix_diag; bp=het_breuschpagan(m.resid,m.model.exog)
    try: reset=linear_reset(m,power=2,use_f=True,cov_type="HC3"); rf,rp=float(reset.fvalue),float(reset.pvalue)
    except Exception: rf,rp=np.nan,np.nan
    diag_rows.append({"Model":name,"N":int(m.nobs),"R2":m.rsquared,"Adjusted_R2":m.rsquared_adj,"Jarque_Bera_p":stats.jarque_bera(m.resid).pvalue,"Breusch_Pagan_p":bp[1],"RESET_F":rf,"RESET_p":rp,"Durbin_Watson":sm.stats.durbin_watson(m.resid),"Max_Cooks_D":cook.max(),"Max_leverage":lev.max(),"Cook_cases_gt_4_over_N":int((cook>4/len(df)).sum())})
    z=model_vif(m); z.insert(0,"Model",name); vif_rows.append(z)
diagnostics=pd.DataFrame(diag_rows); vifs=pd.concat(vif_rows,ignore_index=True)
diagnostics.to_csv(TABLE_DIR/"Table_10_complete_diagnostics.csv",index=False); vifs.to_csv(TABLE_DIR/"Table_S10_model_specific_VIF.csv",index=False)
 
# Supervisor/team clustering audit and clustered-SE sensitivity when an identifier exists.
candidate_exact=["Supervisor_ID","Leader_ID","Team_ID","Dyad_ID","SupervisorID","LeaderID","supervisor_id","leader_id","team_id","dyad_id"]
candidate_regex=[x for x in df.columns if re.search(r"(supervisor|leader|team|dyad).*(id|code)|^(id|code).*(supervisor|leader|team|dyad)",str(x),re.I)]
candidate_all=list(dict.fromkeys(candidate_exact+candidate_regex)); cluster_col=next((x for x in candidate_all if x in df.columns and df[x].nunique(dropna=True)>1),None); cluster_rows=[]
if cluster_col:
    for name,m in principal.items():
        groups=df.loc[m.model.data.row_labels,cluster_col]; cm=m.get_robustcov_results(cov_type="cluster",groups=groups)
        for term,b,se,pv in zip(m.model.exog_names,cm.params,cm.bse,cm.pvalues): cluster_rows.append({"Model":name,"Cluster_variable":cluster_col,"Clusters":groups.nunique(),"Predictor":term,"b":b,"Cluster_SE":se,"p":pv})
else: cluster_rows=[{"Status":"No supervisor, leader, team, or dyad identifier was found; cluster-robust sensitivity could not be estimated."}]
pd.DataFrame(cluster_rows).to_csv(TABLE_DIR/"Table_S13_clustering_audit.csv",index=False)
 
# Participant flow (from recruitment records)
flow=pd.DataFrame({"Stage":["Invited","Usable Time 1","Matched Time 2"],"N":[508,438,391]}); flow.to_csv(TABLE_DIR/"Table_S3_participant_flow.csv",index=False)
# Pearson and Spearman correlation sensitivity.
df[main_vars].corr(method="pearson").to_csv(TABLE_DIR/"Table_S9a_Pearson_correlations.csv")
df[main_vars].corr(method="spearman").to_csv(TABLE_DIR/"Table_S9b_Spearman_correlations.csv")
 
# Age, no-control, and influential-case sensitivity for the main quadratic outcomes.
df["z_Age"]=zscore(df["Age"]); sensitivity=[]
for y,label in outcomes.items():
    for spec,terms in [("Locked controls",REG_CONTROLS+["z_CCB_mean","CCB_sq"]),("Add age",REG_CONTROLS+["z_Age","z_CCB_mean","CCB_sq"]),("No controls",["z_CCB_mean","CCB_sq"])]:
        m,_=fit_hc3(y,terms,spec); sensitivity.append({"Outcome":label,"Specification":spec,"N":int(m.nobs),"Linear_b":m.params.z_CCB_mean,"Linear_p":m.pvalues.z_CCB_mean,"Quadratic_b":m.params.CCB_sq,"Quadratic_p":m.pvalues.CCB_sq,"R2":m.rsquared})
    base=quadratic_models[y]; cooks=base.get_influence().cooks_distance[0]; keep=np.asarray(cooks<=4/len(df))
    design=REG_CONTROLS+["z_CCB_mean","CCB_sq"]; X=sm.add_constant(df.loc[keep,design],has_constant="add")
    robust=sm.OLS(df.loc[keep,y],X).fit(cov_type="HC3")
    sensitivity.append({"Outcome":label,"Specification":"Exclude Cook's D > 4/N","N":int(robust.nobs),"Linear_b":robust.params.z_CCB_mean,"Linear_p":robust.pvalues.z_CCB_mean,"Quadratic_b":robust.params.CCB_sq,"Quadratic_p":robust.pvalues.CCB_sq,"R2":robust.rsquared})
sensitivity=pd.DataFrame(sensitivity); sensitivity.to_csv(TABLE_DIR/"Table_S11_control_sensitivity.csv",index=False)
 
print(f"Wrote {len(list(TABLE_DIR.glob('*.csv')))} tables to {TABLE_DIR}")
