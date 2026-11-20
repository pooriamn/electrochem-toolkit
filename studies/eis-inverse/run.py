"""Reproduce noisy RC/RQ model discrimination and save diagnostic outputs."""
import json
from pathlib import Path
import numpy as np
from electrochem_toolkit import plotstyle
plotstyle.apply()
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from electrochem_toolkit.inverse import Spectrum,response,fit

def main():
    out=Path("study-output/eis");out.mkdir(parents=True,exist_ok=True)
    f=np.logspace(-2,4,60);rng=np.random.default_rng(7)
    truth=response("RQ",f,.12,.8,.025,.82)
    sigma=np.full(len(f),.002)
    z=truth+sigma*(rng.normal(size=len(f))+1j*rng.normal(size=len(f)))
    data=Spectrum(f,z,sigma)
    fits=[fit(data,m) for m in ("RC","RQ")]
    summary={x.model:{"parameters":x.parameters,"standard_errors":x.standard_errors,
                       "weighted_rss":x.weighted_rss,"aicc":x.aicc,
                       "condition_number":x.condition_number,"rank":x.rank,
                       "at_bound":x.at_bound} for x in fits}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summary,indent=2,allow_nan=False))
    fits={x.model:x for x in fits}
    rc,rq=fits["RC"],fits["RQ"]
    style={"RC":dict(color=plotstyle.AMBER,ls="--",lw=1.6),"RQ":dict(color=plotstyle.TEAL,lw=2)}
    obs=dict(s=16,facecolor="white",edgecolor=plotstyle.INK,linewidth=.8,zorder=3)
    plain=FuncFormatter(lambda v,_:f"{v:g}")
    fig,axes=plt.subplots(2,2,figsize=(10.5,7.4),layout="constrained")
    ax=axes[0,0]
    ax.scatter(z.real,-z.imag,label="Synthetic data",**obs)
    for name,x in fits.items(): ax.plot(x.predicted.real,-x.predicted.imag,label=f"{name} fit",**style[name])
    ax.set(xlabel="Re Z / Ω",ylabel="−Im Z / Ω",title="Nyquist");ax.set_aspect("equal");ax.legend(loc="upper right")
    ax=axes[0,1]
    ax.scatter(f,abs(z),**obs)
    for name,x in fits.items(): ax.loglog(f,abs(x.predicted),**style[name])
    ax.set_yticks([.15,.2,.3,.5,.7,1]);ax.yaxis.set_major_formatter(plain);ax.yaxis.set_minor_formatter(FuncFormatter(lambda v,_:""))
    ax.set(xlabel="f / Hz",ylabel="|Z| / Ω",title="Bode magnitude")
    ax=axes[1,0]
    ax.scatter(f,-np.angle(z,deg=True),**obs)
    for name,x in fits.items(): ax.semilogx(f,-np.angle(x.predicted,deg=True),**style[name])
    ax.set(xlabel="f / Hz",ylabel="−phase / °",title="Bode phase")
    ax=axes[1,1];n=len(f)
    ax.axhspan(-2,2,color=plotstyle.TEAL,alpha=.08,lw=0,label="±2σ")
    ax.semilogx(f,rc.standardised_residuals[:n],color=plotstyle.AMBER,lw=1.4,label="RC  real")
    ax.semilogx(f,rc.standardised_residuals[n:],color=plotstyle.AMBER,lw=1.4,ls="--",label="RC  imag")
    ax.semilogx(f,rq.standardised_residuals[:n],color=plotstyle.TEAL,lw=1.4,label="RQ  real")
    ax.semilogx(f,rq.standardised_residuals[n:],color=plotstyle.TEAL,lw=1.4,ls="--",label="RQ  imag")
    ax.axhline(0,color=plotstyle.INK,lw=.6)
    ax.set(xlabel="f / Hz",ylabel="residual / σ",title="Standardised residuals")
    ax.legend(ncol=2,loc="lower left",fontsize=8)
    a,sa=rq.parameters["alpha"],rq.standard_errors["alpha"]
    ax.text(.98,.97,f"ΔAICc (RC − RQ) = {rc.aicc-rq.aicc:,.0f}\nRQ: α = {a:.3f} ± {sa:.3f}  (true 0.82)",
            transform=ax.transAxes,ha="right",va="top",fontsize=8.5,color=plotstyle.INK,linespacing=1.5)
    for ax,l in zip(axes.flat,"abcd"): plotstyle.panel_label(ax,l)
    fig.suptitle("Ideal capacitor or constant-phase element? One noisy spectrum, two fits",x=.01,ha="left",fontsize=13,fontweight="semibold")
    fig.savefig(out/"results.svg",metadata={"Date":None})
    return summary

if __name__=="__main__":main()
