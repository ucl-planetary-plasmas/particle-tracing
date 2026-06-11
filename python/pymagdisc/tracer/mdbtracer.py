from pymagdisc.data.load_data import load_model
from pymagdisc.bfield.magneticField3D import (
    mdiscMagneticField3D,
    dipoleMagneticField3D,
    magneticFieldInterpolator,
)
from pymagdisc.tracer.particle import Particle
from pymagdisc import config
from pymagdisc.tracer.si_constants import eV, c
from pymagdisc.tracer.getbounceperiod import getbounceperiod
from pymagdisc.tracer.getdriftperiod import getdriftperiod
import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from tqdm import tqdm


class MDBTracer:
    """Magnetodisc Boris particle tracer"""

    def __init__(
        self,
        mdfile: str,
        partype: str,
        Ep: float,
        Ri: float,
        ai: float,
        timespec: tuple[float, float, float, float],
        savefile: str = None,
        npertc: int = 5,
    ) -> None:
        """
        Args:
            mdfile (str): magnetodisc mat-file filename.
            partype (str): particle type 'p' -> proton, 'e' -> electron.
            Ep (float): particle energy in MeV.
            Ri (float): initial particle equatorial position (in planet's radius Rp).
            ai (float): initial particle pitch angle (0..180 degrees).
            timespec (list): defined as tmax = sum(timespec*[1,tc,tb,td])
                where 1 is in units of seconds, tc in units of gyroperiod at Ri,
                tb in units of bounce period and td in units of drift period.
            savefile (str, optional): filename to save the simulation data.
                If None, simulation data not saved. Defaults to None.
            npertc (int, optional): number of Boris iterations per gyroperiod. Default to 5.
        """
        self.mdfile = mdfile
        self.partype = partype
        self.Ep = float(Ep)
        self.Ri = float(Ri)
        self.ai = float(ai)
        self.timespec = np.array(timespec, dtype=float)
        self.savefile = savefile
        self.npertc = npertc

    def run_simul(self):
        # load mdsic data
        self.load_data()

        # Create particle
        self.par = Particle(self.partype, self.Ep)

        # get initial conditions
        (
            self.X,
            self.mu,
            self.Tc,
            self.Tb,
            self.Td,
            self.Lm,
            self.rg,
        ) = self.get_init_conditions(R=self.Ri * self.Re, v=self.par.vp, alpha=self.ai)

        # factor mu/gamma^2/m for guiding centre Eq.23 of Ozturk 2012
        self.facdv = self.mu / self.par.gamma**2 / self.par.mp
        self.tspan = [0, sum(self.timespec * np.array([1, self.Tc, self.Tb, self.Td]))]
        print(
            f"tmax={self.tspan[-1]:.10g} s = {self.tspan[-1] / self.Tc:.10g} tc = {self.tspan[-1] / self.Tb:.10g} tb = {self.tspan[-1] / self.Td:.10g} td"
        )

        # start tracing
        results = self.trace(
            self.X,
            self.mdfile,
            self.Ep,
            self.Ri,
            self.ai,
            self.timespec,
            self.Tc,
            self.Tb,
            self.Td,
            self.Lm,
            self.savefile,
            self.npertc,
            self.rg,
        )

        return results

    def load_data(self):
        # load data
        self.MD = load_model(f"{config.PATH_TO_DATA}{self.mdfile}")
        self.Re = float(self.MD["planet"]["r0"])  # Planet's equatorial radius in m
        self.Omega = self.MD["planet"]["Omega"]  # Planet's angular velocity
        self.Omega = self.Omega * interp1d(
            self.MD["disc"]["L"],
            self.MD["disc"]["omega"],
            kind="linear",
            fill_value="extrapolate",
        )(self.Ri)
        print(
            f"Omega planet = {self.MD['planet']['Omega']:.4g} rad/s, "
            f"Omega(Ri) = {self.Omega:.4g} rad/s"
        )
        # Rotation vector in the (x,z)-plane
        # for Saturn 0 deg
        self.tilt = np.deg2rad(0)
        # for Earth and Jupiter 10 deg
        self.tilt = np.deg2rad(10)
        # Rotation vector in the (x,z)-plane
        self.Omp = [self.Omega * np.sin(self.tilt), 0, self.Omega * np.cos(self.tilt)]
        print(f"Dipole tilt = {np.rad2deg(self.tilt):.1f} deg")
        self.B0 = self.MD["planet"]["B0"]  # Magnetic equator field strength in T
        self.B0 = (
            self.B0 * self.MD["v2d"]["Bthdip"][self.MD["dims"]["imu0"] - 1, 0]
        )  # corrected to match Bthdip(mu=0,r=1)
        self.Md = [
            0,
            0,
            self.B0 * self.Re**3,
        ]  # Magnetic moment \mu_0 M/(4\pi) in T m^3
        self.Rm = [0, 0, 0]  # Centered moment

    def get_init_conditions(self, R, v, alpha):
        # initial condition x0,y0,z0,vx0,vy0,vz0
        Xo = [R, 0, 0, 0, v * np.sin(np.deg2rad(alpha)), v * np.cos(np.deg2rad(alpha))]

        #  magnetic field at initial position
        self.interpolator = magneticFieldInterpolator(self.MD)
        B, *_ = mdiscMagneticField3D(
            r=[Xo[0] / self.Re, Xo[1] / self.Re, Xo[2] / self.Re],
            MD=self.MD,
            interpolator=self.interpolator,
        )
        b = np.sqrt(B[3])

        print(f"Ri={R / self.Re:.2f}; pitch angle={alpha:.0f}; B={1e9 * b:.5g} nT")

        # v parallel and perpendicular
        vpar = sum(Xo[3:6] * B[0:3] / b)
        vper = np.sqrt(self.par.V2 - vpar**2)

        # gyro radius at equator
        rg = abs(self.par.gamma / self.par.qOverM * vper / b)
        print(f"rg={rg / self.Re:.4g}Re")

        # first invariant mu
        mu = self.par.gamma**2 * self.par.mp * vper**2 / (2 * b)

        tc, tb, td = self.periods(R, alpha)

        # dipole mirror point latitude
        lm = self.mirrorlat(alpha)

        L = R / self.Re

        # https://farside.ph.utexas.edu/teaching/plasma/lectures1/node22.html
        alphaloss = np.rad2deg(np.arcsin(np.sqrt(1 / np.sqrt(4 * L**6 - 3 * L**5))))
        print(
            f"|alpha| {abs(alpha):.2f}; |pi-alpha| {abs(180 - alpha):.2f}; loss cone {alphaloss:.2f} deg\n"
        )

        if abs(alpha) < alphaloss or abs(180 - alpha) < alphaloss:
            print("WARNING: particle in loss cone!\n")

        return Xo, mu, tc, tb, td, lm, rg

    def trace(
        self, Xo, mdfile, Ep, Ri, ai, timespec, Tc, Tb, Td, Lm, savefile, npertc, rg
    ):
        # Planet's surface
        ts = np.linspace(0, 2 * np.pi, 100)
        Xp = np.cos(ts)
        Yp = np.sin(ts)

        # *dipole* mirror point
        angle = np.deg2rad(Lm.item())  # cleaner than * np.pi / 180
        coslm = np.cos(angle)
        sinlm = np.sin(angle)

        rm = np.array(
            [
                Xo[0] * coslm**3,
                Xo[0] * coslm**2 * sinlm,
                0.0,
            ]
        )

        print("Mirror point = ", rm / self.Re, "Re")

        if np.sqrt(rm[0] ** 2 + rm[1] ** 2) / self.Re < 1:
            print(
                f"Estimated mirror point within planet \n R= {(np.sqrt(rm[0] ** 2 + rm[1] ** 2) / self.Re):.2f} < 1"
            )
            print("MDisc interpolation problem!")

        # B from dipole
        B, gradB = dipoleMagneticField3D(self.Md, self.Rm, rm)
        # Gyro period at dipole mirror point
        omega = (self.par.qOverM / self.par.gamma) * np.sqrt(B[3]).item()
        tm = 2 * np.pi / omega

        # number of iteration such that dt \sim max([tm/3,Tc/5])
        nsteps = int(np.ceil((self.tspan[-1] - self.tspan[0]) / (Tc / npertc)))
        print(self.tspan[0], self.tspan[-1], Tc, nsteps)
        tb = np.linspace(
            self.tspan[0],
            self.tspan[-1],
            nsteps,
        )

        dt = tb[1] - tb[0]

        # tc/tm ratio of gyro period at equator and dipole mirror point
        print(
            f"tm={tm:.4g}, tm/dt={(tm / dt):.4f} tc/tm={(Tc / tm):.4f}, tc/dt={(Tc / dt):.4f}"
        )

        Xb = np.zeros([len(tb), len(Xo)])
        Xf = np.zeros([len(tb), len(Xo)])
        Xgb = np.zeros([len(tb), 3])
        Xgf = np.zeros([len(tb), 3])
        rgb = np.zeros([len(tb), 3])
        rgf = np.zeros([len(tb), 3])

        # tic = time.time()
        Xb[0, :] = self.BorisInit(Xo, dt)  # initial Boris step
        # print("BorisInit time elapsed: ", time.time() - tic)
        Xf[0, :] = Xb[0, :]

        rgb[0, :] = [rg, 0, 0]
        Xgb[0, :] = Xb[0, :3] - rgb[0, :]
        # print(rgb[0,:]/self.Re,Xgb[0,:]/self.Re)

        Xgb[0, :], rgb[0, :] = self.getGyroRadius(
            Xb[0, :3], 0.5 * (Xb[0, 3:6] + Xb[0, 3:6])
        )
        rgb[0, 1:] = 0.0  # initially along x-axis
        # print(rgb[0,:]/self.Re,Xgb[0,:]/self.Re)

        rgf[0, :] = rgb[0, :]
        Xgf[0, :] = Xgb[0, :]

        for i in tqdm(range(len(tb) - 1)):
            Xb[i + 1, :] = self.BorisIter(Xb[i, :], dt)
            Xgb[i + 1, :], rgb[i + 1, :] = self.getGyroRadius(
                Xb[i, :3], 0.5 * (Xb[i, 3:6] + Xb[i + 1, 3:6])
            )

        for i in tqdm(range(len(tb) - 1)):
            Xf[i + 1, :] = self.FullBorisIter(Xf[i, :], dt)
            Xgf[i + 1, :], rgf[i + 1, :] = self.getGyroRadius(
                Xf[i, :3], 0.5 * (Xf[i, 3:6] + Xf[i + 1, 3:6])
            )

        # print(Xgb[:5,:]/self.Re)
        # print(rgb[:5,:]/self.Re)
        # print(Xgf[:5,:]/self.Re)
        # print(rgf[:5,:]/self.Re)

        Zb = Xb[:, 2]
        Rcylb = np.sqrt(Xb[:, 0] ** 2 + Xb[:, 1] ** 2)
        Rtotb = np.sqrt(Xb[:, 0] ** 2 + Xb[:, 1] ** 2 + Xb[:, 2] ** 2)
        Zf = Xf[:, 2]
        Rcylf = np.sqrt(Xf[:, 0] ** 2 + Xf[:, 1] ** 2)
        Rtotf = np.sqrt(Xf[:, 0] ** 2 + Xf[:, 1] ** 2 + Xf[:, 2] ** 2)

        Zgb = Xgb[:, 2]
        Rgcylb = np.sqrt(Xgb[:, 0] ** 2 + Xgb[:, 1] ** 2)
        Rgtotb = np.sqrt(Xgb[:, 0] ** 2 + Xgb[:, 1] ** 2 + Xgb[:, 2] ** 2)
        Zgf = Xf[:, 2]
        Rgcylf = np.sqrt(Xgf[:, 0] ** 2 + Xgf[:, 1] ** 2)
        Rgtotf = np.sqrt(Xgf[:, 0] ** 2 + Xgf[:, 1] ** 2 + Xgf[:, 2] ** 2)

        # Kinetic energy in eV for Full Dynamic Boris
        v2b = np.sum(Xb[:, 3:6] ** 2, axis=1)
        # v2bc = v2b-np.mean(v2b);
        # print(f"v2b {np.mean(v2b):10g} {np.var(v2b, ddof=1):10g}")
        # print(f"v2bc {np.mean(v2bc):10g} {np.var(v2bc, ddof=1):10g}")
        # savemat("v2b_py.mat",{"v2b": v2b, "v2bc": v2bc})
        v2f = np.sum(Xf[:, 3:6] ** 2, axis=1)

        # Eb = 1/2*self.par.mp * v2b / eV;  # classic
        Eb = (1 / np.sqrt(1 - v2b / c**2) - 1) * self.par.mp * c**2 / eV  # relativistic

        Ef1 = (1 / np.sqrt(1 - v2f / c**2) - 1) * self.par.mp * c**2 / eV
        # Only valid for Omega along z-axis
        # Ef2 = -.5*mp*Omega^2*Rcylf.^2/eV; % relativistic
        Ef2 = (
            -0.5
            * self.par.mp
            * (np.dot(self.Omp, self.Omp) * Rtotf**2 - np.dot(Xf[:, :3], self.Omp) ** 2)
            / eV
        )  # relativistic
        # Adjust potential energy of rotation so that minimum is zero
        Ef2 = Ef2 - np.min(Ef2)
        Ef = Ef1 + Ef2

        # mean and variance
        meanEb = np.mean(Eb)
        stdEb = np.std(Eb, ddof=1)
        print(
            f"<Eb> = {(meanEb / 1e6):10g} MeV, std(E) = {stdEb:10g} eV std(E)/<Eb> = {stdEb / meanEb:10g}"
        )
        meanEf = np.mean(Ef)
        stdEf = np.std(Ef, ddof=1)
        print(
            f"<Ef> = {(meanEf / 1e6):10g} MeV, std(E) = {stdEf:10g} eV std(E)/<Ef> = {stdEf / meanEf:10g}"
        )

        # Compute b dot gradB and mu for Full Dynamic Boris
        Bb, gradBb = mdiscMagneticField3D(
            [Xb[:, 0] / self.Re, Xb[:, 1] / self.Re, Xb[:, 2] / self.Re],
            MD=self.MD,
            interpolator=self.interpolator,
            compute_gradB=True,
        )
        bb = np.sqrt(Bb[3])
        BgBb = (Bb[0] * gradBb[0] + Bb[1] * gradBb[1] + Bb[2] * gradBb[2]) / bb
        # Vper2b = V2-(B{1}.*Xb(:,4)+B{2}.*Xb(:,5)+B{3}.*Xb(:,6)).^2./b;
        V2b = Xb[:, 3] ** 2 + Xb[:, 4] ** 2 + Xb[:, 5] ** 2
        Vparb = (Bb[0] * Xb[:, 3] + Bb[1] * Xb[:, 4] + Bb[2] * Xb[:, 5]) / bb
        Vperb = np.sqrt(V2b - Vparb**2)
        Vper2b = V2b - Vparb**2
        # Instantaneous first invariant Eq.14
        muib = (self.par.gamma**2) * self.par.mp * Vper2b / (2 * bb)
        # pitch angle
        aib = 180 - np.rad2deg(np.arctan2(Vperb, Vparb))

        Bf, gradBf = mdiscMagneticField3D(
            [Xf[:, 0] / self.Re, Xf[:, 1] / self.Re, Xf[:, 2] / self.Re],
            MD=self.MD,
            interpolator=self.interpolator,
            compute_gradB=True,
        )
        bf = np.sqrt(Bf[3])
        BgBf = (Bf[0] * gradBf[0] + Bf[1] * gradBf[1] + Bf[2] * gradBf[2]) / bf
        # Vper2b = V2-(B{1}.*Xb(:,4)+B{2}.*Xb(:,5)+B{3}.*Xb(:,6)).^2./b;
        V2f = Xf[:, 3] ** 2 + Xf[:, 4] ** 2 + Xf[:, 5] ** 2
        Vparf = (Bf[0] * Xf[:, 3] + Bf[1] * Xf[:, 4] + Bf[2] * Xf[:, 5]) / bf
        Vperf = np.sqrt(V2f - Vparf**2)
        Vper2f = V2f - Vparf**2
        # Instantaneous first invariant Eq.14
        muif = (self.par.gamma**2) * self.par.mp * Vper2f / (2 * bf)
        # pitch angle
        aif = 180 - np.rad2deg(np.arctan2(Vperf, Vparf))

        # Compute latitude and longitude
        latb = np.rad2deg(np.arctan2(Zb, Rcylb))  # atan(Rcyl/Z)
        lonb = np.rad2deg(np.unwrap(np.arctan2(Xb[:, 1], Xb[:, 0])))  # atan(X/Y)
        latf = np.rad2deg(np.arctan2(Zf, Rcylf))  # atan(Rcyl/Z)
        lonf = np.rad2deg(np.unwrap(np.arctan2(Xf[:, 1], Xf[:, 0])))  # atan(X/Y)

        # L-Shell with L
        t = np.pi * np.linspace(-0.5, 0.5, 100)
        # estimate L from instantaneous L
        # coslat = Rcyl./R; Ls = R./coslat.^2/Re;
        Lsb = Rtotb**3 / Rcylb**2 / self.Re
        Le = np.mean(Lsb)
        print(f"L estimated = {Le:.2f}\n")

        r = Le * np.cos(t) ** 2
        ir = np.argwhere(r >= 1).reshape(-1)
        xe = r[ir] * np.cos(t[ir])
        ye = r[ir] * np.sin(t[ir])

        print("*** Bounce period Boris")
        # identify zero crossings
        izc = np.argwhere(latb[0:-2] * latb[1:-1] < 0).reshape(-1)
        Tbe = 2 * np.mean(np.diff(tb[izc]))
        # initial period value for fit
        # average of dipole Tb and Tbe with weight 1 and 2
        Tbi = (Tb + 2 * Tbe) / 3
        # exact period
        Tbi = Tbe
        print(f"Init Tbd={Tb:.2f} Tbe={Tbe:.2f} Tbi={Tbi:.2f} ****")
        tbb, dtbb, lmb, fitlatb = getbounceperiod(tb, latb, 2 * np.pi / Tbi)

        print("*** Bounce period FullBoris")
        # identify zero crossings
        izc = np.argwhere(latf[0:-2] * latf[1:-1] < 0).reshape(-1)
        Tbe = 2 * np.mean(np.diff(tb[izc]))
        # initial period value for fit
        # average of dipole Tb and Tbe with weight 1 and 2
        Tbi = (Tb + 2 * Tbe) / 3
        # exact period
        Tbi = Tbe
        print(f"Init Tbd={Tb:.2f} Tbe={Tbe:.2f} Tbi={Tbi:.2f} ****")
        tbf, dtbf, lmf, fitlatf = getbounceperiod(tb, latf, 2 * np.pi / Tbi)

        print("*** Drift period Boris")
        tdb, dtdb, fitlonb = getdriftperiod(tb, lonb, 2 * np.pi / tbb, 360)

        print("*** Drift period FullBoris")
        tdf, dtdf, fitlonf = getdriftperiod(tb, lonf, 2 * np.pi / tbf, 360)

        return {
            "mdfile": self.mdfile,
            "Re": self.Re,
            "Be": self.Be,
            "partype": self.partype,
            "par": self.par,
            "Ep": Ep,
            "Ri": Ri,
            "ai": ai,
            "Omega": self.Omega,
            "Omp": self.Omp,
            "timespec": timespec,
            "npertc": self.npertc,
            "Tc": Tc,
            "Tb": Tb,
            "Td": Td,
            "Lm": Lm,
            "rg": rg,
            "tb": tb,
            "Xb": Xb,
            "Xf": Xf,
            "Xgb": Xgb,
            "rgb": rgb,
            "Xgf": Xgf,
            "rgf": rgf,
            "Zb": Zb,
            "Rcylb": Rcylb,
            "Rtotb": Rtotb,
            "Eb": Eb,
            "muib": muib,
            "aib": aib,
            "latb": latb,
            "lonb": lonb,
            "Bb": Bb,
            "bb": bb,
            "gradBb": gradBb,
            "Vparb": Vparb,
            "Vperb": Vperb,
            "Vper2b": Vper2b,
            "Zf": Zf,
            "Rcylf": Rcylf,
            "Rtotf": Rtotf,
            "Ef": Ef,
            "muif": muif,
            "aif": aif,
            "latf": latf,
            "lonf": lonf,
            "Bf": Bf,
            "bf": bf,
            "gradBf": gradBf,
            "Vparf": Vparf,
            "Vperf": Vperf,
            "Vper2f": Vper2f,
            "tbb": tbb,
            "dtbb": dtbb,
            "lmb": lmb,
            "fitlatb": fitlatb,
            "tdb": tdb,
            "dtdb": dtdb,
            "fitlonb": fitlonb,
            "tbf": tbf,
            "dtbf": dtbf,
            "lmf": lmf,
            "fitlatf": fitlatf,
            "tdf": tdf,
            "dtdf": dtdf,
            "fitlonf": fitlonf,
        }

    def FullBorisIter(self, r, h):
        r = np.array(r)

        x, y, z, vx, vy, vz = r

        B, _ = mdiscMagneticField3D(
            [x / self.Re, y / self.Re, z / self.Re],
            MD=self.MD,
            interpolator=self.interpolator,
        )

        Bxyz = np.array([B[0], B[1], B[2]])

        fac = self.par.qOverM / self.par.gamma

        # T = -[Bx,By,Bz]'/sqrt(B{4}*tan(\theta/2.0) = qB/m\Delta{t}/2;
        T = fac * Bxyz * h / 2
        S = 2.0 * T / (1 + np.dot(T, T))

        V = np.array([vx, vy, vz])
        R = np.array([x, y, z])

        # https://en.wikipedia.org/wiki/Centrifugal_force
        # Centrifugal force dv/dt = -\omega x (\omega x r)
        # Only valid for Omega along z-axis
        # dV = Omega^2*[x,y,0]'*h/2;
        dV = -np.cross(self.Omp, np.cross(self.Omp, R)) * h / 2
        # https://en.wikipedia.org/wiki/Coriolis_force
        # Coriolis force dv/dt = - 2 \omega x dr/dt
        # Only valid for Omega along z-axis
        # dV = -2*Omega*[-vy,vx,0]'*h/2;
        # dV = -2*cross(Omp,V)*h/2;
        # Centrifugal and Coriolis forces
        # Only valid for Omega along z-axis
        # dV = (Omega*[x,y,0]'*h/2 -2*Omega*[-vy,vx,0]')*h/2;
        # dV = -cross(Omp,cross(Omp,R))*h/2-2*cross(Omp,V)*h/2;

        vm = V + dV
        vp = vm + np.cross(vm, T)
        vp = vm + np.cross(vp, S)
        V = vp + dV

        R = R + V * h

        rn = np.concatenate([R, V])

        return rn

    def BorisIter(self, r, h):
        r = np.array(r)

        x, y, z, vx, vy, vz = r

        B, _ = mdiscMagneticField3D(
            [x / self.Re, y / self.Re, z / self.Re],
            MD=self.MD,
            interpolator=self.interpolator,
        )

        Bxyz = np.array([B[0], B[1], B[2]])

        fac = self.par.qOverM / self.par.gamma

        # T = -[Bx,By,Bz]'/sqrt(B{4}*tan(\theta/2.0) = qB/m\Delta{t}/2;
        T = fac * Bxyz * h / 2
        S = 2.0 * T / (1 + np.dot(T, T))

        V = np.array([vx, vy, vz])
        R = np.array([x, y, z])

        v = V + np.cross(V, T)
        V = V + np.cross(v, S)

        R = R + V * h

        rn = np.concatenate([R, V])

        return rn

    def BorisInit(self, r, h):
        # https://physics.stackexchange.com/questions/296863/how-to-initialize-bootstrap-the-boris-algorithm

        # dt is now -dt/2 (to go backward in time for dt/2)
        h = -h / 2

        r = np.array(r)

        x, y, z, vx, vy, vz = r

        B, _ = mdiscMagneticField3D(
            [x / self.Re, y / self.Re, z / self.Re],
            MD=self.MD,
            interpolator=self.interpolator,
        )

        Bxyz = np.array([B[0], B[1], B[2]])

        fac = self.par.qOverM / self.par.gamma

        # T = -[Bx,By,Bz]'/sqrt(B{4}*tan(\theta/2.0) = qB/m\Delta{t}/2;
        T = fac * Bxyz * h / 2
        S = 2.0 * T / (1 + np.dot(T, T))

        V = np.array([vx, vy, vz])
        R = np.array([x, y, z])

        v = V + np.cross(V, T)
        V = V + np.cross(v, S)

        rn = np.concatenate([R, V])

        return rn

    def periods(self, R, alpha):
        self.B, *_ = mdiscMagneticField3D(
            r=[R / self.Re, 0, 0], MD=self.MD, interpolator=self.interpolator
        )
        self.Be, *_ = mdiscMagneticField3D(
            r=[self.Re / self.Re, 0, 0], MD=self.MD, interpolator=self.interpolator
        )

        self.B = np.sqrt(self.B[3])
        self.Be = np.sqrt(self.Be[3])

        # Angular gyro frequency
        Omegac = abs(self.par.qOverM) / self.par.gamma * self.B

        # Gyro period
        tc = 2 * np.pi / Omegac

        # Bouncing period
        # factor Jupiter 1.2478 = Re/c*sqrt(2)*3.7
        # factor Saturn  1.0521 = Re/c*sqrt(2)*3.7
        fac = self.Re / c * np.sqrt(2) * 3.7
        tb = (
            fac
            * (R / self.Re)
            * c
            / self.par.vp
            * (1 - 0.4635 * np.sin(np.deg2rad(alpha)) ** 0.75)
        )
        # https://farside.ph.utexas.edu/teaching/plasma/lectures1/node22.html
        # tb = R*sqrt(2)/vp*(3.7-1.6*sind(alpha));

        # Drift period
        td = (
            2
            * np.pi
            * abs(self.par.qOverM)
            * self.Be
            * self.Re**3
            / self.par.vp**2
            / R
            * (1 - 1 / 3 * (np.sin(np.deg2rad(alpha))) ** 0.62)
        )
        # https://farside.ph.utexas.edu/teaching/plasma/lectures1/node23.html
        # td = pi*qp*Be*Re^3/(3*EMeV*R)/(0.35+0.15*sind(alpha));
        # td = 1.05/(EMeV/(1e6*eV))/(R/Re)/(1+0.43*sind(alpha))*3600;

        assert td > tb > tc, (
            "Drift period should be largerthan bounce, and bounce larger than gyro."
        )
        print(f"tc={tc:.2f} s, tb={tb:.2f} s, td={td:.2f} s\n")

        return tc, tb, td

    def mirrorlat(self, alpha):
        # numerically

        def myfun(x, a):
            x = np.deg2rad(x[0])  # extract scalar from fsolve input
            a = np.deg2rad(a)
            return [np.cos(x) ** 6 - np.sin(a) ** 2 * np.sqrt(1 + 3 * np.sin(x) ** 2)]

        lm = fsolve(func=myfun, x0=[45], args=(alpha,))
        lm_val = lm[0]

        print(f"pitcheq={alpha:.2f}; lm={lm_val:.2f} deg\n")

        return abs(lm_val)

    def getGyroRadius(self, Rp, Vp):
        tol = 1e-10

        Rp = np.asarray(Rp, dtype=float).reshape(3)
        Vp = np.asarray(Vp, dtype=float).reshape(3)

        rg2 = 0.0
        relerr = np.inf

        # initial guiding centre position
        Xg = Rp.copy()

        fac = -self.par.gamma / self.par.qOverM

        while relerr > tol:
            x, y, z = Xg

            # magnetic field at guiding centre position
            B, _ = mdiscMagneticField3D(
                [x / self.Re, y / self.Re, z / self.Re],
                MD=self.MD,
                interpolator=self.interpolator,
            )

            Bxyz = np.array([B[0], B[1], B[2]])
            b2 = B[3]

            # parallel velocity projection
            Vppar = (np.dot(Bxyz, Vp) / b2) * Bxyz

            # perpendicular velocity
            Vpper = Vp - Vppar

            # gyro radius
            rg = fac * np.cross(Vpper, Bxyz) / b2

            oldrg2 = rg2
            rg2 = np.dot(rg, rg)

            relerr = abs(rg2 - oldrg2) / max(rg2, np.finfo(float).eps)

            Xg = Rp - rg

        return Xg, rg
