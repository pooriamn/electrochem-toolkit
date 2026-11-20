"""Reproduce a dynamic stack calculation and export energy-balance diagnostics."""
import csv
import json
from pathlib import Path
from electrochem_toolkit import plotstyle
plotstyle.apply()
import matplotlib.pyplot as plt
from electrochem_toolkit.electrothermal import Stack,simulate

def main():
    out=Path("study-output/electrothermal");out.mkdir(parents=True,exist_ok=True)
    rows=simulate(Stack(),[(300,.1,0),(900,.6,0),(900,1.0,0),(600,.2,0),(600,0,0)],sample_s=10)
    with (out/"timeseries.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    summary={"peak_temperature_K":max(r["temperature_K"] for r in rows),
             "final_temperature_K":rows[-1]["temperature_K"],
             "max_abs_balance_residual_J":max(abs(r["balance_residual_J"]) for r in rows),
             "final_dc_energy_J":rows[-1]["dc_energy_J"],
             "final_chemical_energy_J":rows[-1]["chemical_energy_J"],
             "final_heat_loss_J":rows[-1]["heat_loss_J"],
             "final_thermal_storage_J":rows[-1]["thermal_storage_J"]}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summary,indent=2,allow_nan=False))
    t=[r["time_s"]/60 for r in rows]
    col=lambda k:[r[k] for r in rows]
    fig,axes=plt.subplots(2,2,figsize=(10.5,7.4),layout="constrained",sharex=True)
    ax=axes[0,0];ax.step(t,col("current_density_A_cm2"),where="post",color=plotstyle.SLATE)
    ax.fill_between(t,col("current_density_A_cm2"),step="post",color=plotstyle.SLATE,alpha=.08)
    ax.set(ylabel="j / A cm⁻²",title="Load programme")
    ax=axes[0,1];temp=[x-273.15 for x in col("temperature_K")];ax.plot(t,temp,color=plotstyle.TEAL)
    k=max(range(len(temp)),key=temp.__getitem__)
    ax.annotate(f"peak {temp[k]:.1f} °C",(t[k],temp[k]),xytext=(8,-4),textcoords="offset points",fontsize=8.5,va="top")
    ax.set(ylabel="temperature / °C",title="Stack temperature")
    ax=axes[1,0];v,etn=col("cell_voltage_V"),col("E_tn_V")
    ax.fill_between(t,etn,v,where=[a>b for a,b in zip(v,etn)],color=plotstyle.AMBER,alpha=.15,lw=0,label="net heat release")
    ax.plot(t,v,color=plotstyle.INK,label="cell")
    ax.plot(t,etn,color=plotstyle.AMBER,ls="--",lw=1.4,label="thermoneutral")
    ax.plot(t,col("E_rev_V"),color=plotstyle.TEAL,lw=1.4,label="reversible")
    ax.set(xlabel="time / min",ylabel="voltage / V",title="Cell voltage against its thermodynamic limits");ax.legend(loc="upper right",fontsize=8)
    ax=axes[1,1]
    for field,label,c in [("dc_energy_J","DC input",plotstyle.INK),("chemical_energy_J","stored as H₂ (ΔH)",plotstyle.TEAL),
                          ("heat_loss_J","rejected heat",plotstyle.AMBER),("thermal_storage_J","held in stack",plotstyle.SLATE)]:
        ax.plot(t,[x/1e6 for x in col(field)],label=label,color=c)
    res=max(abs(x) for x in col("balance_residual_J"))
    ax.text(.03,.6,f"max balance residual {res:.1e} J",transform=ax.transAxes,ha="left",fontsize=8.5,color=plotstyle.NICKEL)
    ax.set(xlabel="time / min",ylabel="energy / MJ",title="Cumulative energy balance");ax.legend(loc="upper left",fontsize=8)
    for ax,l in zip(axes.flat,"abcd"): plotstyle.panel_label(ax,l)
    fig.suptitle("Illustrative 20-cell electrolyser stack under a step load: voltage, heat and energy closure",x=.01,ha="left",fontsize=13,fontweight="semibold")
    fig.savefig(out/"results.svg",metadata={"Date":None})
    return summary

if __name__=="__main__":main()
