#!/usr/bin/env python3
# =============================================================================
# refugee_wsbm.py — refugee roles from a weighted stochastic blockmodel
#
#   python analysis/refugee_wsbm.py [output_dir]           (default: output)
#   Needs graph-tool (conda-forge: graph-tool-base), numpy, pandas.
#   Reads <output_dir>/unhcr_wsbm/edges.csv (written by analysis/refugee_wsbm_edges.R)
#   and <output_dir>/unhcr/roles_by_country_period.csv (the REGE roles).
#
# Why a blockmodel. Regular equivalence groups countries whose ties go to
# equivalent others. A stochastic blockmodel (SBM) groups countries whose ties,
# and here also the sizes of their flows, follow the same pattern with respect
# to the groups. It is the probabilistic form of blockmodelling (Holland,
# Laskey & Leinhardt 1983). With REGE on weighted refugee flows, the result
# depends on when the iterations are stopped (analysis/refugee_rege_weighted.R).
# The SBM has no such choice: the number of roles, the roles themselves and the
# model variant are all chosen by one criterion, the description length (DL,
# Peixoto 2017), i.e. the number of bits needed to describe the network with
# the model. The fit that compresses the data most is preferred; more roles are
# kept only if they save more bits than they cost.
#
# Network: one node per country-year with any refugee flow (the stacked
# network of the pipeline, 1990-2020); one directed edge per flow of at least
# one refugee, weighted by the number of refugees. Country-years of different
# years are not linked, but they can share a role if they are linked to the
# roles in the same way, so roles are comparable over time.
#
# Models compared by DL (each from 20 independent fits, best kept):
#   binary              edges only, no flow sizes
#   real-normal (log)   log(flow) normally distributed within each pair of roles
#   discrete-geometric  flow counts geometric within each pair of roles
# each with and without degree correction (degree correction lets countries of
# one role differ in their number of partners), flat and nested (hierarchical).
# DLs of weighted models are comparable with each other on the scale of the
# flows (the log-normal DL is converted with its Jacobian, sum of log flows).
#
# For the chosen model, the posterior distribution of partitions is sampled by
# MCMC (merge-split moves), and the roles are the posterior mode
# (PartitionModeState); each country-year also gets the probability of its role.
#
# Robustness: (1) agreement between independent fits; (2) ties of at least 10 or
# 100 refugees; (3) each year fitted separately; (4) the other weight models;
# compared with the adjusted Rand index (ARI; 1 = same roles, 0 = chance).
#
# Writes to <output_dir>/unhcr_wsbm/:
#   model_selection.csv      DL of every model variant (bits), number of roles
#   roles_wsbm.csv           role of every country-year, its posterior probability
#   role_profiles.csv        size, sending/receiving, flows, REGE roles in each role
#   role_image.csv           edges and mean log flow between roles (image graph)
#   robustness.csv           ARI between the main roles and each alternative
#   by_year.csv              each year fitted separately: roles, ARI with main
# =============================================================================
import sys, os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import graph_tool.all as gt

OUT = sys.argv[1] if len(sys.argv) > 1 else "output"
D = os.path.join(OUT, "unhcr_wsbm")
SEED = 2024
N_FITS = int(os.environ.get("WSBM_FITS", 20))
gt.seed_rng(SEED); np.random.seed(SEED)

edges = pd.read_csv(os.path.join(D, "edges.csv"))
rege = pd.read_csv(os.path.join(OUT, "unhcr", "roles_by_country_period.csv"))

def ari(a, b):
    """adjusted Rand index between two labelings of the same items"""
    a = pd.factorize(pd.Series(a))[0]; b = pd.factorize(pd.Series(b))[0]
    t = pd.crosstab(a, b).values; c2 = lambda x: x * (x - 1) / 2
    s = c2(t).sum(); sa = c2(t.sum(1)).sum(); sb = c2(t.sum(0)).sum(); n = c2(t.sum())
    e = sa * sb / n; m = (sa + sb) / 2
    return float((s - e) / (m - e)) if m != e else 1.0

def build(e):
    """directed graph of country-years with flow and log-flow edge properties"""
    nodes = sorted(set(e.origin + "_" + e.period.astype(str)) | set(e.destination + "_" + e.period.astype(str)))
    idx = {n: i for i, n in enumerate(nodes)}
    g = gt.Graph(directed=True); g.add_vertex(len(nodes))
    src = [idx[o + "_" + str(p)] for o, p in zip(e.origin, e.period)]
    dst = [idx[d + "_" + str(p)] for d, p in zip(e.destination, e.period)]
    g.add_edge_list(np.column_stack([src, dst]))
    g.ep.flow = g.new_ep("int", vals=np.round(e.flow.values).astype(int))
    g.ep.logw = g.new_ep("double", vals=np.log(e.flow.values))
    return g, nodes

def state_args(g, w, dc):
    """arguments of the base state: degree correction and, for weighted models, the
    edge covariate. Weighted models need WeightedBlockState (in this graph-tool
    version unknown keyword arguments of BlockState are silently ignored)."""
    a = dict(deg_corr=dc)
    if w == "real-normal (log)": a.update(rec=[g.ep.logw], rec_types=["real-normal"])
    if w == "discrete-geometric": a.update(rec=[g.ep.flow], rec_types=["discrete-geometric"])
    return a

def fit(g, w, dc, nested, n=N_FITS):
    """best of n independent minimisations, refined by greedy merge-split sweeps"""
    base = gt.BlockState if w == "binary" else gt.WeightedBlockState
    best = None
    for _ in range(n):
        if nested:
            s = gt.minimize_nested_blockmodel_dl(g, base_state=base, base_state_args=state_args(g, w, dc))
        else:
            s = gt.minimize_blockmodel_dl(g, state=base, state_args=state_args(g, w, dc))
        for _ in range(10): s.multiflip_mcmc_sweep(beta=np.inf, niter=10)
        if best is None or s.entropy() < best.entropy(): best = s
    return best

def labels(s, nested):
    st = s.get_levels()[0] if nested else s
    return np.array(st.get_blocks().a)

def n_roles(s, nested):
    return len(np.unique(labels(s, nested)))

LN2 = np.log(2)
g, nodes = build(edges)
jac = float(g.ep.logw.a.sum())           # nats: converts the log-normal DL to the scale of the flows
print(f"Network: {g.num_vertices()} country-years, {g.num_edges()} flows", flush=True)

# ---- model selection ------------------------------------------------------------
rows, fits = [], {}
for w in ["binary", "real-normal (log)", "discrete-geometric"]:
    for dc in [True, False]:
        for nested in [False, True]:
            s = fit(g, w, dc, nested)
            dl = s.entropy()                                   # nats
            dl_flows = dl + (jac if w == "real-normal (log)" else 0)
            rows.append(dict(weights=w, degree_corrected=dc, nested=nested, roles=n_roles(s, nested),
                             levels=(len([l for l in s.get_levels() if len(np.unique(l.b.a)) > 1]) if nested else 1),
                             dl_bits=dl / LN2, dl_bits_flow_scale=dl_flows / LN2))
            fits[(w, dc, nested)] = s
            print(rows[-1], flush=True)
ms = pd.DataFrame(rows)
wt = ms[ms.weights != "binary"]
ms["chosen"] = False
best = wt.dl_bits_flow_scale.idxmin(); ms.loc[best, "chosen"] = True
ms.to_csv(os.path.join(D, "model_selection.csv"), index=False)
W, DC, NEST = ms.loc[best, ["weights", "degree_corrected", "nested"]]
print("Chosen:", W, "degree-corrected" if DC else "", "nested" if NEST else "flat", flush=True)

# ---- posterior of the chosen model ---------------------------------------------
s = fits[(W, DC, NEST)]
bs = []
def collect(st):
    bs.append(st.get_bs() if NEST else st.b.a.copy())
gt.mcmc_equilibrate(s, force_niter=int(os.environ.get("WSBM_NITER", 2000)), mcmc_args=dict(niter=10), callback=collect)
pm = gt.PartitionModeState(bs, nested=NEST, converge=True)
if NEST:
    bmax = pm.get_max_nested(); lab = np.array(bmax[0])
    upper = [np.array(b) for b in bmax[1:]]
else:
    lab = np.array(pm.get_max(g).a); upper = []
marg = pm.get_marginal(g)
prob = np.array([marg[v].a[lab[int(v)]] / marg[v].a.sum() if marg[v].a.sum() > 0 else np.nan for v in g.vertices()])
# relabel roles 0..K-1 by size
order = pd.Series(lab).value_counts().index.tolist(); lab = np.array([order.index(x) + 1 for x in lab])

df = pd.DataFrame({"node": nodes, "role": lab, "posterior_prob": prob})
df[["iso3", "period"]] = df.node.str.rsplit("_", n=1, expand=True); df["period"] = df.period.astype(int)
# coarser levels of the hierarchy (nested model), mapped down to the nodes
if NEST:
    cur = np.array(bmax[0])
    for lv, b in enumerate(bmax[1:], start=1):
        b = np.array(b); cur = b[cur]
        if len(np.unique(cur)) <= 1: break
        o = pd.Series(cur).value_counts().index.tolist(); df[f"role_level{lv}"] = [o.index(x) + 1 for x in cur]
df = df.merge(rege[["iso3", "period", "role"]].rename(columns={"role": "rege_role"}), on=["iso3", "period"], how="left")
st = g.degree_property_map("out").a, g.degree_property_map("in").a
fo = np.zeros(g.num_vertices()); fi = np.zeros(g.num_vertices())
for e_, f in zip(g.get_edges(), g.ep.flow.a): fo[e_[0]] += f; fi[e_[1]] += f
df["out_partners"], df["in_partners"], df["out_flow"], df["in_flow"] = st[0], st[1], fo, fi
df.to_csv(os.path.join(D, "roles_wsbm.csv"), index=False)
print("Roles:", df.role.value_counts().sort_index().to_dict(), "| mean posterior prob", round(np.nanmean(prob), 3), flush=True)

def profile(df, col):
    return df.groupby(col).apply(lambda x: pd.Series(dict(
        n=len(x), sends=(x.out_partners > 0).mean(), receives=(x.in_partners > 0).mean(),
        median_out_flow=x.out_flow.median(), median_in_flow=x.in_flow.median(),
        median_partners=(x.out_partners + x.in_partners).median(),
        rege_sender_only=(x.rege_role == 1).mean(), rege_two_way=(x.rege_role == 2).mean(), rege_receiver_only=(x.rege_role == 3).mean(),
        posterior_prob=x.posterior_prob.mean(),
        examples=" ".join(x.assign(tot=x.out_flow + x.in_flow).sort_values("tot", ascending=False).iso3.drop_duplicates().head(10))))).reset_index()
prof = profile(df, "role"); prof.to_csv(os.path.join(D, "role_profiles.csv"), index=False)
for c in [c for c in df.columns if c.startswith("role_level")]:
    profile(df, c).to_csv(os.path.join(D, f"role_profiles_{c}.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(prof.round(2).to_string(), flush=True)

# image graph: edges and mean log flow between roles
E = pd.DataFrame(g.get_edges(), columns=["s", "t"]); E["from"] = lab[E.s]; E["to"] = lab[E.t]; E["logw"] = g.ep.logw.a
img = E.groupby(["from", "to"]).agg(edges=("logw", "size"), mean_log_flow=("logw", "mean")).reset_index()
img.to_csv(os.path.join(D, "role_image.csv"), index=False)

# ---- robustness ----------------------------------------------------------------
rob = []
main = pd.Series(lab, index=nodes)
rob.append(dict(check="REGE roles of the pipeline (binary, 3 iterations)", ari=ari(df.role, df.rege_role), roles=3))
for lv in [c for c in df.columns if c.startswith("role_level")]:
    rob.append(dict(check=f"REGE roles vs blockmodel {lv}", ari=ari(df[lv], df.rege_role), roles=df[lv].nunique()))
# independent fits of the chosen model
for k in range(10):
    s2 = fit(g, W, DC, NEST, n=1)
    rob.append(dict(check="independent fit of the chosen model", ari=ari(lab, labels(s2, NEST)), roles=n_roles(s2, NEST)))
# the other model variants
for key, s2 in fits.items():
    if key == (W, DC, NEST): continue
    rob.append(dict(check=f"model: {key[0]}, {'degree-corrected' if key[1] else 'not degree-corrected'}, {'nested' if key[2] else 'flat'}",
                    ari=ari(lab, labels(s2, key[2])), roles=n_roles(s2, key[2])))
# tie thresholds
for thr in [10, 100]:
    e2 = edges[edges.flow >= thr]; g2, n2 = build(e2); s2 = fit(g2, W, DC, NEST, n=5)
    l2 = pd.Series(labels(s2, NEST), index=n2); common = l2.index.intersection(main.index)
    rob.append(dict(check=f"ties of at least {thr} refugees", ari=ari(main[common], l2[common]), roles=n_roles(s2, NEST), nodes=len(common)))
rob = pd.DataFrame(rob); rob.to_csv(os.path.join(D, "robustness.csv"), index=False)
print(rob.round(3).to_string(), flush=True)

# each year separately
yr = []
for y, e2 in edges.groupby("period"):
    g2, n2 = build(e2); s2 = fit(g2, W, DC, NEST, n=5)
    l2 = pd.Series(labels(s2, NEST), index=n2)
    yr.append(dict(year=y, nodes=len(n2), edges=len(e2), roles_year_fit=n_roles(s2, NEST),
                   roles_main_in_year=main[n2].nunique(), ari_with_main=ari(main[n2], l2),
                   ari_rege_year=ari(df.set_index("node").loc[n2, "rege_role"], l2)))
yr = pd.DataFrame(yr); yr.to_csv(os.path.join(D, "by_year.csv"), index=False)
print(yr.round(3).to_string(), flush=True)
print("Done:", D)
