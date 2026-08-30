#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HumidClimatologyEngine v10 analytical comparison toolkit.

Reads v10 period NetCDF files and produces reproducible CSV/PNG/HTML summaries.
Designed for post-processing only: it never modifies source NetCDF files.
"""
from __future__ import annotations
import argparse, json, math, hashlib
from pathlib import Path
from datetime import datetime, timezone
import numpy as np

ENGINE_VERSION = "10.0.0"
PERIODS = ("1981_1990","1991_2000","2001_2010","2011_2020","1981_2020")
VARIABLES = ("rh","e","r","q")
LEVELS = {"L1":0,"L2":1,"L3":9}
LEVEL_RANGES = {"L1":(0,1),"L2":(1,9),"L3":(9,33)}

def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def pct_change(a,b):
    den=np.abs(b)
    out=np.full_like(a,np.nan,dtype=float)
    mask=np.isfinite(a)&np.isfinite(b)&(den>0)
    out[mask]=100*(a[mask]-b[mask])/den[mask]
    return out

def find_period_file(root:Path, period:str)->Path:
    candidates=list(root.rglob(f"*{period}*.nc"))
    # Prefer main statistics file over diagnostics/bivariate/histogram when ambiguous.
    candidates=[p for p in candidates if 'diagnostic' not in p.name.lower() and 'bivariate' not in p.name.lower() and 'hist' not in p.name.lower()]
    if not candidates: raise FileNotFoundError(f"No main NetCDF found for period {period} under {root}")
    candidates.sort(key=lambda p:(len(p.name),p.name))
    return candidates[0]

def load_dataset(path:Path):
    try:
        import xarray as xr
    except Exception as e:
        raise RuntimeError("xarray is required for the analysis toolkit") from e
    return xr.open_dataset(path, engine='netcdf4', decode_times=True, mask_and_scale=True, cache=False)

def load_core(root:Path,period:str,var:str,level:str):
    p=find_period_file(root,period)
    with load_dataset(p) as ds:
        name=f"mean_{var}"
        if name not in ds: raise KeyError(f"{name} not found in {p}")
        da=ds[name].isel(level_bin=LEVELS[level])
        n=ds[f"n_{var}"].isel(level_bin=LEVELS[level])
        return np.asarray(da.values,dtype=float), np.asarray(n.values,dtype=float), np.asarray(ds.latitude.values), np.asarray(ds.longitude.values), p

def weighted_domain_summary(values,n):
    mask=np.isfinite(values)&np.isfinite(n)&(n>0)
    if not np.any(mask): return {"mean":math.nan,"min":math.nan,"max":math.nan,"cells":0,"support":0}
    w=n[mask]; x=values[mask]
    return {"mean":float(np.average(x,weights=w)),"min":float(np.nanmin(x)),"max":float(np.nanmax(x)),"cells":int(mask.sum()),"support":int(w.sum())}

def csv_write(path,rows,headers):
    import csv
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=headers); w.writeheader(); w.writerows(rows)

def make_plots(outdir, summary_rows, maps, l1_change):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return []
    out=[]
    # Decadal overview plot per variable
    for var in VARIABLES:
        xs=[r['period'] for r in summary_rows if r['variable']==var and r['level']=='L1']
        ys=[float(r['weighted_mean']) for r in summary_rows if r['variable']==var and r['level']=='L1']
        if len(xs):
            fig=plt.figure(figsize=(9,5)); ax=fig.add_subplot(111); ax.plot(xs,ys,marker='o'); ax.set_title(f"{var.upper()} L1 weighted domain mean - v10"); ax.set_xlabel('Period'); ax.set_ylabel(var); ax.grid(alpha=.25); fig.tight_layout(); p=outdir/f"decadal_{var}_L1_mean.png"; fig.savefig(p,dpi=180); plt.close(fig); out.append(str(p))
    # Change maps for 2011-2020 vs 1981-1990, one per variable
    for var,(arr,lat,lon) in l1_change.items():
        fig=plt.figure(figsize=(9,5)); ax=fig.add_subplot(111); im=ax.imshow(arr,origin='upper',aspect='auto'); fig.colorbar(im,ax=ax,label=f"% change {var}"); ax.set_title(f"L1 {var.upper()} percent change: 2011-2020 vs 1981-1990"); ax.set_xlabel('Longitude index'); ax.set_ylabel('Latitude index'); fig.tight_layout(); p=outdir/f"change_{var}_2011_2020_vs_1981_1990_L1.png"; fig.savefig(p,dpi=180); plt.close(fig); out.append(str(p))
    return out

def run(args):
    root=Path(args.input_root); out=Path(args.output_dir); out.mkdir(parents=True,exist_ok=True)
    rows=[]; maps={}; meta={"engine_version":ENGINE_VERSION,"created_utc":datetime.now(timezone.utc).isoformat(),"input_root":str(root)}
    for level in args.levels:
        for var in VARIABLES:
            for period in PERIODS:
                values,n,lat,lon,p=load_core(root,period,var,level)
                s=weighted_domain_summary(values,n)
                rows.append({"period":period,"variable":var,"level":level,"weighted_mean":s['mean'],"min":s['min'],"max":s['max'],"cells":s['cells'],"support":s['support'],"source":str(p)})
                meta.setdefault('sources',{})[period]=str(p)
    headers=list(rows[0].keys())
    csv_write(out/'domain_summary.csv',rows,headers)
    # Build primary change maps
    for var in VARIABLES:
        a,_,_,_,_=load_core(root,'1981_1990',var,'L1'); b,_,lat,lon,_=load_core(root,'2011_2020',var,'L1'); maps[var]=pct_change(b,a)
    # compact difference tables
    change_rows=[]
    for var in VARIABLES:
        for level in args.levels:
            lookup={(r['period'],r['variable'],r['level']):r for r in rows if r['variable']==var and r['level']==level}
            base=lookup[('1981_1990',var,level)]; latest=lookup[('2011_2020',var,level)]; full=lookup[('1981_2020',var,level)]
            for metric in ('weighted_mean','min','max'):
                bv=float(base[metric]); lv=float(latest[metric]); fv=float(full[metric])
                change_rows.append({"variable":var,"level":level,"metric":metric,"1981_1990":bv,"2011_2020":lv,"absolute_change":lv-bv,"percent_change":(100*(lv-bv)/abs(bv) if np.isfinite(bv) and bv!=0 else math.nan),"FULL":fv})
    csv_write(out/'decadal_change_2011_2020_vs_1981_1990.csv',change_rows,list(change_rows[0].keys()))
    plots=make_plots(out,rows,maps,{k:(v,np.asarray(lat),np.asarray(lon)) for k,v in maps.items()})
    meta['plots']=plots; meta['sha256_sources']={period:sha256_file(find_period_file(root,period)) for period in PERIODS}
    (out/'analysis_manifest.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    # Simple HTML index
    html=['<!doctype html><html><head><meta charset="utf-8"><title>HumidClimatologyEngine v10 Analysis</title></head><body>',f'<h1>HumidClimatologyEngine v{ENGINE_VERSION} Analysis</h1>',f'<p>Created {meta["created_utc"]}</p>','<h2>Decadal change</h2>','<table border="1"><tr><th>Variable</th><th>Level</th><th>Metric</th><th>1981-90</th><th>2011-20</th><th>Absolute</th><th>Percent</th></tr>']
    for r in change_rows:
        html.append('<tr>'+''.join(f'<td>{r[k]}</td>' for k in ('variable','level','metric','1981_1990','2011_2020','absolute_change','percent_change'))+'</tr>')
    html.append('</table><h2>Figures</h2>')
    for p in plots:
        rel=Path(p).name; html.append(f'<p><img src="{rel}" style="max-width:900px"></p>')
    html.append('</body></html>'); (out/'analysis_report.html').write_text('\n'.join(html),encoding='utf-8')
    print(json.dumps({"status":"PASS","output":str(out),"plots":plots},indent=2))

def main():
    ap=argparse.ArgumentParser(description='HumidClimatologyEngine v10 analytical comparison toolkit')
    ap.add_argument('--input-root',required=True); ap.add_argument('--output-dir',required=True)
    ap.add_argument('--levels',nargs='+',choices=('L1','L2','L3'),default=['L1','L2','L3'])
    args=ap.parse_args(); run(args)

if __name__=='__main__': main()
