from pathlib import Path

from scipy.io import savemat

from pymagdisc.tracer.mdbtracer import MDBTracer

if __name__ == "__main__":
    mdfile = "jup_mdisc_kh3e7_rmp90.mat"  # magnetodisc mat-file
    partype = "p"                         # proton
    Ep = 100                              # energy 100 MeV
    Ri = 4                                # initial equatorial distance 4*Rp
    ai = 30                               # initial pitch angle 30 degrees
    timespec = [0, 0, 2, 0]               # run for 2 dipole bounce periods
    savefile = "my_jup_mdisc_sim"         # name for result mat-file
    npertc = 5                            # number of Boris iterations per gyroperiod

    # Create particle tracer object
    tracer = MDBTracer(
        mdfile=mdfile,
        partype=partype,
        Ep=Ep,
        Ri=Ri,
        ai=ai,
        timespec=timespec,
        savefile=savefile,
        npertc=npertc,
    )

    # Run simulation
    results = tracer.run_simul()

    # Create output directory if needed
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)

    # Save output data
    output_path = results_dir / f"{savefile}.mat"
    savemat(output_path, results)

    print(f"Results saved to: {output_path.relative_to(Path('.'))}")
