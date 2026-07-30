
"""
Independent numerical mirror used to validate the MATLAB equations in an
environment where MATLAB/Octave is unavailable. It is not required to run
the MATLAB project.
"""
from pathlib import Path
import json, math, csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
e = 1.602176634e-19
hbar = 1.054571817e-34
kB = 1.380649e-23
eps0 = 8.8541878128e-12
m0 = 9.1093837015e-31
meV = e*1e-3
nm = 1e-9

def softplus(x):
    return np.maximum(x, 0) + np.log1p(np.exp(-np.abs(x)))

def trap_weights(x):
    w = np.zeros_like(x)
    w[0] = (x[1]-x[0])/2
    w[-1] = (x[-1]-x[-2])/2
    w[1:-1] = (x[2:]-x[:-2])/2
    return w

def solve_sp(Nz=700, B=10.0, T=300.0):
    m = 0.067*m0
    Lz = 5.0*nm
    domain = 15.0*nm
    z = np.linspace(-domain/2, domain/2, Nz)
    dz = z[1]-z[0]
    u = z/Lz
    Vconf = 228.0*meV*u**2*(0.30*u**6-1)
    Nd_sheet = 1e17
    width = 2*nm
    mask = np.abs(z) <= width/2
    Nd_z = np.zeros_like(z)
    Nd_z[mask] = Nd_sheet/np.trapezoid(mask.astype(float), z)
    VH = np.zeros_like(z)
    EF_old = np.nan
    converged = False
    for it in range(220):
        Veff = Vconf + VH + e**2*B**2*z**2/(2*m)
        a = hbar**2/(2*m*dz**2)
        H = diags([-a*np.ones(Nz-3), 2*a+Veff[1:-1],
                   -a*np.ones(Nz-3)], [-1,0,1], format='csr')
        E, P = eigsh(H, k=4, which='SA')
        idx = np.argsort(E)
        E, P = E[idx], P[:,idx]
        P /= np.sqrt(np.sum(P**2, axis=0)*dz)
        Psi = np.zeros((Nz,4))
        Psi[1:-1,:] = P
        pref = m*kB*T/(np.pi*hbar**2)
        def neutrality(EF):
            return pref*np.sum(softplus((EF-E)/(kB*T))) - Nd_sheet
        EF = brentq(neutrality, E.min()-100*kB*T,
                    E.max()+100*kB*T, xtol=1e-30)
        Ni = pref*softplus((EF-E)/(kB*T))
        n_z = (np.abs(Psi)**2)@Ni
        rhs = e**2/(13.18*eps0)*(Nd_z-n_z)
        D2 = diags([np.ones(Nz-3), -2*np.ones(Nz-2),
                    np.ones(Nz-3)], [-1,0,1], format='csr')
        VH_new = np.zeros_like(z)
        VH_new[1:-1] = spsolve(D2, dz**2*rhs[1:-1])
        VH_mix = 0.2*VH_new + 0.8*VH
        perr = np.max(np.abs(VH_mix-VH))/max(np.max(np.abs(VH_mix)),1e-3*meV)
        eferr = np.inf if not np.isfinite(EF_old) else abs(EF-EF_old)/meV
        VH = VH_mix
        if perr <= 1e-4 and eferr <= 1e-6:
            converged = True
            break
        EF_old = EF

    # final solve
    Veff = Vconf + VH + e**2*B**2*z**2/(2*m)
    a = hbar**2/(2*m*dz**2)
    H = diags([-a*np.ones(Nz-3), 2*a+Veff[1:-1],
               -a*np.ones(Nz-3)], [-1,0,1], format='csr')
    E, P = eigsh(H,k=4,which='SA')
    idx=np.argsort(E); E,P=E[idx],P[:,idx]
    P/=np.sqrt(np.sum(P**2,axis=0)*dz)
    Psi=np.zeros((Nz,4)); Psi[1:-1]=P
    pref=m*kB*T/(np.pi*hbar**2)
    EF=brentq(lambda x: pref*np.sum(softplus((x-E)/(kB*T)))-Nd_sheet,
              E.min()-100*kB*T,E.max()+100*kB*T,xtol=1e-30)
    Ni=pref*softplus((EF-E)/(kB*T))
    return dict(z=z,Vconf=Vconf,VH=VH,Veff=Veff,E=E,Psi=Psi,EF=EF,
                Ni=Ni,Nd_sheet=Nd_sheet,converged=converged,
                iterations=it+1,perr=perr,eferr=eferr,T=T,B=B)

def form_factor(sp,qz):
    product=np.conj(sp["Psi"][:,0])*sp["Psi"][:,1]
    return np.trapezoid(product[:,None]*np.exp(1j*np.outer(sp["z"],qz)),
                        sp["z"],axis=0)

def binned_lorentzian(Egrid, centers, weights, order, gamma):
    dE=Egrid[1]-Egrid[0]
    density=np.zeros_like(Egrid)
    idx=np.rint((centers.ravel()-Egrid[0])/dE).astype(int)
    val=weights.ravel()/(order*dE)
    ok=(idx>=0)&(idx<Egrid.size)&np.isfinite(val)&(val!=0)
    np.add.at(density,idx[ok],val[ok])
    off=(np.arange(Egrid.size)-Egrid.size//2)*dE
    g=gamma/order
    kernel=g/np.pi/(off**2+g**2)
    return np.convolve(density,kernel,mode="same")*dE

def compute_spectrum(sp,Nq=60,dE_meV=0.15):
    m=0.067*m0
    Egrid=np.arange(2.0,160.0+dE_meV/2,dE_meV)*meV
    qz=np.linspace(0,2.5/nm,Nq)
    qp=np.linspace(0,1.5/nm,Nq)
    wz,wp=trap_weights(qz),trap_weights(qp)
    QP,QZ=np.meshgrid(qp,qz)
    WP,WZ=np.meshgrid(wp,wz)
    q=np.sqrt(QP**2+QZ**2)
    lB=np.sqrt(hbar/(e*sp["B"]))
    J2=np.exp(-0.5*lB**2*QP**2)
    measure=2/(2*np.pi)**2*QP*WP*WZ*J2
    I2=np.abs(form_factor(sp,qz))**2
    I2grid=I2[:,None]
    ne3d=1.0e18*1.0e6
    qd=np.sqrt(ne3d*e**2/(13.18*eps0*kB*sp["T"]))
    dE=sp["E"][1]-sp["E"][0]
    dip=np.trapezoid(np.conj(sp["Psi"][:,0])*sp["z"]*sp["Psi"][:,1],sp["z"])
    weights_mb=np.exp(-(sp["E"]-sp["E"].min())/(kB*sp["T"]))
    pop=sp["Nd_sheet"]*weights_mb/weights_mb.sum()
    E0=4.5e5
    a0=7.5*nm
    pref=E0**2*np.sqrt(13.18)/(8*np.pi)*(2*np.pi/hbar)
    Mrad=e**2*E0**2*abs(dip)**2/4
    results={"optical":{},"piezoelectric":{}}
    total=np.zeros_like(Egrid)

    for mech in results:
        if mech=="optical":
            hw=np.full_like(q,36.25*meV)
            C2V=e**2*(36.25*meV)/(2*eps0)*(1/10.89-1/13.18)*q**2/(q**2+qd**2)**2
            gamma=0.8*meV
        else:
            s=5.22e3
            hw=hbar*s*q
            C2V=0.006*hbar*e**2*s/(2*13.18*eps0)*q**3/(q**2+qd**2)**2
            gamma=0.6*meV
        N=np.zeros_like(q)
        mask=hw>0
        N[mask]=1/np.expm1(hw[mask]/(kB*sp["T"]))
        base=measure*C2V*I2grid
        mtot=np.zeros_like(Egrid)
        for ell in [1,2,3]:
            dressing=QP**(2*ell)/(2**(2*ell)*math.factorial(ell)**2)
            em=binned_lorentzian(Egrid,(dE+hw)/ell,
                                 base*dressing*(N+1),ell,gamma)
            ab=binned_lorentzian(Egrid,(dE-hw)/ell,
                                 base*dressing*N,ell,gamma)
            scale=pref*pop[0]*Mrad*a0**(2*ell)/Egrid**2
            results[mech][ell]={"emission":scale*em,
                                "absorption":scale*ab,
                                "total":scale*(em+ab)}
            mtot+=results[mech][ell]["total"]
        results[mech]["total"]=mtot
        total+=mtot
    return Egrid,total,results,dict(qd=qd,dip=dip,deltaE=dE)

def peak_metrics(x,y):
    i=int(np.argmax(y))
    peak=float(y[i])
    half=peak/2
    left=np.where(y[:i+1]<=half)[0]
    right=np.where(y[i:]<=half)[0]
    if left.size and right.size and right[0]>0:
        il=left[-1]
        ir=i+right[0]
        xl=np.interp(half,[y[il],y[il+1]],[x[il],x[il+1]])
        xr=np.interp(half,[y[ir-1],y[ir]],[x[ir-1],x[ir]])
        fwhm=float(xr-xl)
    else:
        fwhm=float("nan")
    return {"peak_x_meV":float(x[i]/meV),"fwhm_meV":fwhm/meV,
            "hwhm_meV":fwhm/(2*meV)}

def run():
    rows=[]
    mirrors={}
    for name,Nz,Nq,step in [("quick",500,38,0.25),("standard",900,70,0.12)]:
        sp=solve_sp(Nz=Nz)
        Eg,total,res,meta=compute_spectrum(sp,Nq=Nq,dE_meV=step)
        mirrors[name]=(sp,Eg,total,res,meta)
        rows.append({
            "profile":name,
            "converged":sp["converged"],
            "iterations":sp["iterations"],
            "E0_meV":sp["E"][0]/meV,
            "E1_meV":sp["E"][1]/meV,
            "DeltaE01_meV":(sp["E"][1]-sp["E"][0])/meV,
            "EF_meV":sp["EF"]/meV,
            "charge_rel_error":abs(sp["Ni"].sum()-sp["Nd_sheet"])/sp["Nd_sheet"],
            "norm_max_error":float(np.max(np.abs(np.trapezoid(abs(sp["Psi"])**2,sp["z"],axis=0)-1))),
            "total_peak_meV":Eg[np.argmax(total)]/meV
        })

    sp,Egrid,total,res,meta=mirrors["standard"]
    report={
        "runtime_note":"Independent Python mirror; MATLAB/Octave unavailable in build container.",
        "profiles":rows,
        "standard":{
            "energies_meV":(sp["E"]/meV).tolist(),
            "EF_meV":float(sp["EF"]/meV),
            "deltaE01_meV":float(meta["deltaE"]/meV),
            "dipole01_nm":float(meta["dip"]/nm),
            "qd_inv_nm":float(meta["qd"]*nm),
            "total_metrics":peak_metrics(Egrid,total),
            "components":{}
        }
    }
    for mech in ["optical","piezoelectric"]:
        report["standard"]["components"][mech]={}
        for ell in [1,2,3]:
            report["standard"]["components"][mech][str(ell)]={
                "emission":peak_metrics(Egrid,res[mech][ell]["emission"]),
                "absorption":peak_metrics(Egrid,res[mech][ell]["absorption"])
            }

    (ROOT/"validation_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    with (ROOT/"convergence_profiles.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys())
        w.writeheader(); w.writerows(rows)

    plt.figure(figsize=(7,5))
    plt.plot(sp["z"]/nm,sp["Vconf"]/meV,linestyle="--",label="UA")
    plt.plot(sp["z"]/nm,sp["Veff"]/meV,label="Ueff")
    scale=40
    for i in range(4):
        rho=abs(sp["Psi"][:,i])**2
        plt.plot(sp["z"]/nm,sp["E"][i]/meV+scale*rho/rho.max(),label=f"E{i}")
    plt.xlabel("z (nm)"); plt.ylabel("Energy (meV)")
    plt.title("Independent validation: electronic structure")
    plt.grid(True); plt.legend()
    plt.tight_layout(); plt.savefig(ROOT/"structure_validation.png",dpi=180); plt.close()

    scale=total.max()
    plt.figure(figsize=(7,5))
    plt.plot(Egrid/meV,total/scale,label="Total")
    plt.plot(Egrid/meV,res["optical"]["total"]/scale,label="Optical phonon")
    plt.plot(Egrid/meV,res["piezoelectric"]["total"]/scale,label="Piezoelectric phonon")
    plt.xlabel("Photon energy (meV)"); plt.ylabel("Normalized MOAP")
    plt.title("Independent validation: total spectrum")
    plt.grid(True); plt.legend()
    plt.tight_layout(); plt.savefig(ROOT/"spectrum_validation.png",dpi=180); plt.close()

    plt.figure(figsize=(8,5))
    for ell in [1,2,3]:
        plt.plot(Egrid/meV,res["optical"][ell]["total"]/scale,label=f"Optical {ell}PA")
    plt.xlabel("Photon energy (meV)"); plt.ylabel("Normalized MOAP")
    plt.title("Independent validation: optical photon orders")
    plt.grid(True); plt.legend()
    plt.tight_layout(); plt.savefig(ROOT/"orders_validation.png",dpi=180); plt.close()

    print(json.dumps(report,indent=2))

if __name__=="__main__":
    run()
