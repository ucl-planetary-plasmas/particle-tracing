from pymagdisc.data.load_data import load_model
from pymagdisc.bfield.magneticField3D import (
    mdiscMagneticField3D,
    dipoleMagneticField3D,
    magneticFieldInterpolator,
)
from pymagdisc.tracer.particle import Particle
from pymagdisc import config
from pymagdisc.tracer.si_constants import e, eV, c, mp, me
from pymagdisc.tracer.getbounceperiod import getbounceperiod
from pymagdisc.tracer.getdriftperiod import getdriftperiod
import numpy as np
from scipy.optimize import fsolve
from tqdm import tqdm
import time


class MDBTracer:
    """Magnetodisc Boris particle tracer"""

    def __init__(
        self,
        mdfile: str,
        partype: str,
        Ep: float,
        Ri: float,
        ai: float,
        timespec: list,
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
        self.Ep = Ep
        self.Ri = Ri
        self.ai = ai
        self.timespec = timespec
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
            self.rlarm,
        ) = self.get_init_conditions(R=self.Ri * self.Re, v=self.par.vp, alpha=self.ai)

        # factor mu/gamma^2/m for guiding centre Eq.23 of Ozturk 2012
        self.facdv = self.mu / self.par.gamma ** 2 / self.par.mp
        self.tspan = [0, sum(self.timespec * np.array([1, self.Tc, self.Tb, self.Td]))]
        print(
            f"tmax={self.tspan[-1]:.2f} s = {self.tspan[-1]/self.Tc:.2f} tc = {self.tspan[-1]/self.Tb:.2f} tb = {self.tspan[-1]/self.Td:.2f} td"
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
            self.rlarm,
        )

        return results

    def load_data(self):
        # load data
        self.MD = load_model(f"{config.PATH_TO_DATA}{self.mdfile}")
        self.Re = self.MD["planet"]["r0"]  # Planet's equatorial radius in m
        self.B0 = self.MD["planet"]["B0"]  # Magnetic equator field strength in T
        self.B0 = (
            self.B0 * self.MD["v2d"]["Bthdip"][self.MD["dims"]["imu0"] - 1, 0]
        )  # corrected to match Bthdip(mu=0,r=1)
        self.Md = [
            0,
            0,
            self.B0 * self.Re ** 3,
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

        print(f"R={R/self.Re:.2f}; pitch angle={alpha:.0f}; B={1e9*b:.2f} nT\n")

        # v parallel and perpendicular
        vpar = sum(Xo[3:6] * B[0:3] / b)
        vper = np.sqrt(self.par.V2 - vpar ** 2)

        # gyro radius
        rho = abs(self.par.gamma / self.par.qOverM * vper / b)
        rlarm = self.par.rho

        # first invariant mu
        mu = self.par.gamma ** 2 * mp * vper ** 2 / (2 * b)

        tc, tb, td = self.periods(R, alpha)

        # dipole mirror point latitude
        lm = self.mirrorlat(alpha)

        L = R / self.Re

        # https://farside.ph.utexas.edu/teaching/plasma/lectures1/node22.html
        alphaloss = (
            np.arcsin(np.sqrt(1 / np.sqrt(4 * L ** 6 - 3 * L ** 5))) * 180 / np.pi
        )
        print(
            f"|alpha| {abs(alpha):.2f}; |pi-alpha| {abs(180-alpha):.2f}; loss cone {alphaloss:.2f} deg\n"
        )

        if abs(alpha) < alphaloss or abs(180 - alpha) < alphaloss:
            print("WARNING: particle in loss cone!\n")

        return Xo, mu, tc, tb, td, lm, rlarm

    def trace(
        self, Xo, mdfile, Ep, Ri, ai, timespec, Tc, Tb, Td, Lm, savefile, npertc, rlarm
    ):

        # Planet's surface
        ts = np.linspace(0, 2 * np.pi, 100)
        Xp = np.cos(ts)
        Yp = np.sin(ts)

        # *dipole* mirror point
        rm = [
            float(Xo[0] * np.cos(Lm * np.pi / 180) ** 3),
            float(Xo[0] * (np.cos(Lm * np.pi / 180) ** 2) * np.sin(Lm * np.pi / 180)),
            float(0),
        ]

        print("Mirror point = ", rm)

        if np.sqrt(rm[0] ** 2 + rm[1] ** 2) / self.Re < 1:
            print(
                f"Estimated mirror point within planet \n R= {(np.sqrt(rm[0]**2+rm[1]**2)/self.Re):.2f} < 1"
            )
            print("MDisc interpolation problem!")

        # B from dipole
        B, gradB = dipoleMagneticField3D(self.Md, self.Rm, rm)
        # Gyro period at dipole mirror point
        tm = float(2 * np.pi / (self.par.qOverM / self.par.gamma * np.sqrt(B[3])))

        # number of iteration such that dt \sim max([tm/3,Tc/5])
        tb = np.linspace(
            self.tspan[0],
            self.tspan[-1],
            int(np.ceil(np.diff(self.tspan) / [Tc / npertc])),
        )

        dt = float(tb[1] - tb[0])

        # tc/tm ratio of gyro period at equator and dipole mirror point
        print(
            f"tm={tm:.4f}, tm/dt={(tm/dt):.4f} tc/tm={(Tc/tm):.4f}, tc/dt={(Tc/dt):.4f}"
        )

        Xb = np.zeros([len(tb), len(Xo)])

        tic = time.time()
        Xb[0, :] = self.BorisIter(Xo, dt, initStep=True)  # initial Boris step
        print("BorisInit time elapsed: ", time.time() - tic)

        for i in tqdm(range(len(tb) - 1)):
            Xb[i + 1, :] = self.BorisIter(Xb[i, :], dt)

        Zb = Xb[:, 2]
        Rcylb = np.sqrt(Xb[:, 0] ** 2 + Xb[:, 1] ** 2)
        Rtotb = np.sqrt(Xb[:, 0] ** 2 + Xb[:, 1] ** 2 + Xb[:, 2] ** 2)

        # Kinetic energy in eV for Full Dynamic Boris
        v2b = np.sum(Xb[:, 3:6] ** 2, axis=1)

        # Eb = 1/2*mp*v2b/eV;  % classic
        Eb = (1 / np.sqrt(1 - v2b / c ** 2) - 1) * mp * c ** 2 / eV  # relativistic

        # mean and variance
        meanEb = np.mean(Eb)
        stdEb = np.std(Eb)
        print(f"<Eb> = {(meanEb/1e6):.2f} MeV, std(E) = {stdEb:.2e} eV\n")

        # Compute b dot gradB and mu for Full Dynamic Boris
        B, gradB = mdiscMagneticField3D(
            [Xb[:, 0] / self.Re, Xb[:, 1] / self.Re, Xb[:, 2] / self.Re],
            MD=self.MD,
            interpolator=self.interpolator,
            compute_gradB=True,
        )
        b = np.sqrt(B[3])
        BgBb = (B[0] * gradB[0] + B[1] * gradB[1] + B[2] * gradB[2]) / b
        # Vper2b = V2-(B{1}.*Xb(:,4)+B{2}.*Xb(:,5)+B{3}.*Xb(:,6)).^2./b;
        V2b = Xb[:, 3] ** 2 + Xb[:, 4] ** 2 + Xb[:, 5] ** 2
        Vparb = (B[0] * Xb[:, 3] + B[1] * Xb[:, 4] + B[2] * Xb[:, 5]) / b
        Vperb = np.sqrt(V2b - Vparb ** 2)
        Vper2b = self.par.V2 - Vparb ** 2
        # Instantaneous first invariant Eq.14
        muib = (self.par.gamma ** 2) * mp * Vper2b / (2 * b)

        # Compute latitude and longitude
        latb = np.arctan2(Zb, Rcylb) * 180 / np.pi  # atan(Rcyl/Z)
        lonb = np.arctan2(Xb[:, 1], Xb[:, 0]) * 180 / np.pi  # atan(X/Y);

        pitchan = 180 - np.arctan2(Vperb, Vparb) * 180 / np.pi

        # L-Shell with L
        t = np.linspace(-np.pi / 2, np.pi / 2, 100)
        # estimate L from instantaneous L
        # coslat = Rcyl./R; Ls = R./coslat.^2/Re;
        Lsb = Rtotb ** 3 / Rcylb ** 2 / self.Re
        Le = np.mean(Lsb)
        print(f"L estimated = {Le:.2f}\n")

        r = Le * np.cos(t) ** 2
        ir = np.argwhere(r >= 1).reshape(-1)
        xe = r[ir] * np.cos(t[ir])
        ye = r[ir] * np.sin(t[ir])

        # identify zero crossings
        izc = np.argwhere(latb[0:-2] * latb[1:-1] < 0).reshape(-1)
        Tbe = 2 * np.mean(np.diff(tb[izc]))
        # Taking average of Tb and Tbe with weight 1 and 2
        Tbi = (Tb + 2 * Tbe) / 3
        Tbi = Tbe
        print(f"**** Tbd={Tb:.2f} Tbe={Tbe:.2f} Tbi={Tbi:.2f} ****\n")
        tbb, dtbb, lmb, fitlat = getbounceperiod(tb, latb, 2 * np.pi / Tbi)
        tdb, dtdb, fitlon = getdriftperiod(tb, lonb, 2 * np.pi / tbb, 360)

        return {
            "mdfile": self.mdfile,
            "Re": self.Re,
            "Be": self.Be,
            "Ep": Ep,
            "Ri": Ri,
            "ai": ai,
            "timespec": timespec,
            "Tc": Tc,
            "Tb": Tb,
            "Td": Td,
            "Lm": Lm,
            "tb": tb,
            "Xb": Xb,
            "Zb": Zb,
            "Rcylb": Rcylb,
            "Rtotb": Rtotb,
            "Eb": Eb,
            "muib": muib,
            "latb": latb,
            "lonb": lonb,
            "tbb": tbb,
            "dtbb": dtbb,
            "lmb": lmb,
            "fitlat": fitlat,
            "tdb": tdb,
            "dtdb": dtdb,
            "fitlon": fitlon,
            "pitchan": pitchan,
            "rho": self.par.rho,
        }

    def BorisIter(self, r, h, initStep=False):

        if initStep:
            # dt is now -dt/2 (to go backward in time for dt/2)
            h = -h / 2

        r = np.array(r)
        rn = np.zeros(r.shape)

        [x, y, z, vx, vy, vz] = [r[0], r[1], r[2], r[3], r[4], r[5]]

        B, _ = mdiscMagneticField3D(
            [x / self.Re, y / self.Re, z / self.Re],
            MD=self.MD,
            interpolator=self.interpolator,
        )

        Bxyz = np.array([B[0], B[1], B[2]])

        fac = self.par.qOverM / self.par.gamma

        # T = -[Bx,By,Bz]'/sqrt(B{4}*tan(\theta/2.0) = qB/m\Delta{t}/2;
        T = fac * Bxyz * h / 2
        Ttr = np.transpose(T)
        S = 2.0 * T / (1 + Ttr @ T)

        V = np.array([vx, vy, vz])
        R = np.array([x, y, z])

        v = V + np.cross(V, T)
        V = V + np.cross(v, S)

        if not initStep:
            R = R + V * h

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
            * self.Re ** 3
            / self.par.vp ** 2
            / R
            * (1 - 1 / 3 * (np.sin(np.deg2rad(alpha))) ** 0.62)
        )
        # https://farside.ph.utexas.edu/teaching/plasma/lectures1/node23.html
        # td = pi*qp*Be*Re^3/(3*EMeV*R)/(0.35+0.15*sind(alpha));
        # td = 1.05/(EMeV/(1e6*eV))/(R/Re)/(1+0.43*sind(alpha))*3600;

        assert (
            td > tb > tc
        ), "Drift period should be largerthan bounce, and bounce larger than gyro."
        print(f"tc={tc:.2f} s, tb={tb:.2f} s, td={td:.2f} s\n")

        return tc, tb, td

    def mirrorlat(self, alpha):

        # numerically

        myfun = lambda x, a: np.cos(np.deg2rad(x)) ** 6 - np.sin(
            np.deg2rad(a)
        ) ** 2 * np.sqrt(1 + 3 * np.sin(np.deg2rad(x)) ** 2)
        lm = fsolve(func=myfun, x0=[45], args=(alpha))
        print(f"pitcheq={alpha:.2f}; lm={lm[0]:.2f} deg\n")

        return abs(lm)
