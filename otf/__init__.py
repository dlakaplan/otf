try:
    from presto.prepfold import pfd
except ImportError:
    print("PRESTO python libraries not available")
import numpy as np
from astropy import units as u, constants as c
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, Longitude, Latitude
from scipy.optimize import curve_fit
from astropy.io import fits
import scipy.ndimage
import os
import matplotlib.pyplot as plt


def extract_pointing(fitsfile, suffix="pointing", format="ecsv"):
    """
    outfile = extract_pointing(fitsfile, suffix="pointing", format="ecsv")

    extracts the OFFS_SUB, RA_SUB, DEC_SUB columns from a PSRFITS file and writes to a separate astropy.table
    format defines the extension
    """
    with fits.open(fitsfile) as f:
        t = Table(
            [
                Column(f[1].data["OFFS_SUB"], name="OFFS_SUB"),
                Column(f[1].data["RA_SUB"], name="RA_SUB"),
                Column(f[1].data["DEC_SUB"], name="DEC_SUB"),
            ]
        )
        outfile = os.path.splitext(fitsfile)[0] + f".{suffix}." + format
        t.write(outfile, overwrite=True)
        return outfile


class PulseProfile:
    def __init__(self, *args, **kwargs):
        """
        PulseProfile(*args, **kwargs)

        usage:
        PulseProfile(gaussianfitfile)

        or:
        PulseProfile(fwhm=fwhm, pha=pha, [ampl=ampl], [const=const])

        """
        if len(args) == 1:
            # load via file
            self.load(args[0])
        else:
            if not "fwhm" in kwargs:
                self.fwhm = np.array([])
            else:
                self.fwhm = np.atleast_1d(kwargs["fwhm"])
            if not "pha" in kwargs:
                self.pha = np.array([])
            else:
                self.pha = np.atleast_1d(kwargs["pha"])
            if not "ampl" in kwargs:
                self.ampl = np.ones_like(self.fwhm)
            else:
                self.ampl = np.atleast_1d(kwargs["ampl"])
            if not "const" in kwargs:
                self.const = 0
            else:
                self.const = kwargs["const"]
            assert len(self.pha) == len(self.ampl) == len(self.fwhm), (
                "Number of phases (%d), amplitudes (%d), and FWHMs (%d) are not the same"
                % (len(self.pha), len(self.ampl), len(self.fwhm))
            )

    @classmethod
    def fit_gaussian_profile(cls, x, y, fwhm=0.1, pha=None, ampl=None):
        """
        fit_gaussian_profile(x, y, fwhm=0.1, pha=None, ampl=None)
        returns a PulseProfile object resulting from a single Gaussian fit to y(x)
        if not provided, initial guesses for pha and ampl are from where y==y.max()

        """
        if pha is None:
            pha = x[y == y.max()].mean()
        if ampl is None:
            ampl = y.max()
        popt, _ = curve_fit(gaussian, x, y, p0=[fwhm, pha, ampl])
        return cls(fwhm=popt[0], pha=popt[1], ampl=popt[2])

    def __len__(self):
        return len(self.fwhm)

    @property
    def len(self):
        """
        returns self.len as a property
        """
        return len(self.fwhm)

    @property
    def sigma(self):
        """
        returns sigma (converted from FWHM)
        """
        return self.fwhm / (2 * np.sqrt(2 * np.log(2)))

    def __call__(self, x, normalize=True, ignoreconst=True):
        """
        y = object(x, normalize=True, ignoreconst=True)

        return the PulseProfile evaluated at x
        if normalize is False then ignores the amplitudes
        """
        y = np.zeros_like(x)
        for i in range(len(self.fwhm)):
            if normalize:
                y += gaussian(x, self.fwhm[i], self.pha[i], self.ampl[i])
            else:
                y += gaussian(x, self.fwhm[i], self.pha[i], 1)
        if not ignoreconst:
            y += self.const
        return y

    def load(self, gaussians_file):
        """
        load(gaussians_file)
        loads in the Gaussian fit file and stores results

        """
        with open(gaussians_file) as f:
            lines = f.readlines()
        phass = []
        ampls = []
        fwhms = []
        for line in lines:
            if line.lstrip().startswith("const"):
                const = float(line.split()[2])
            if line.lstrip().startswith("phas"):
                phass.append(float(line.split()[2]))
            if line.lstrip().startswith("ampl"):
                ampls.append(float(line.split()[2]))
            if line.lstrip().startswith("fwhm"):
                fwhms.append(float(line.split()[2]))

        assert len(phass) == len(ampls) == len(fwhms), (
            "Number of phases (%d), amplitudes (%d), and FWHMs (%d) are not the same in '%s'!"
            % (len(phass), len(ampls), len(fwhms), gaussians_file)
        )

        self.pha = np.asarray(phass)
        self.ampl = np.asarray(ampls)
        self.fwhm = np.asarray(fwhms)
        self.const = const


def gaussian(x, fwhm, x0, ampl):
    """
    y = gaussian(x, fwhm, x0, ampl)
    returns a Gaussian with the specified parameters
    """
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return ampl * np.exp(-((x - x0) ** 2) / 2 / sigma ** 2)


class OTF_Scan:
    def __init__(
        self, pfdfile, pointingfile, profilefile=None, offpulse=None, autoroll=True
    ):
        """
        OTF_Scan(pfdfile, pointingfile, profilefile=None, offpulse=None, autoroll=True)

        create an OTF_Scan object
        needs a pfd file and a pointing file (a Table which has columns OFFS_SUB, RA_SUB, DEC_SUB)

        if specified will also read in a file with the results of pygaussfit to determine the profile

        offpulse will be pairs [min,max] where min and max are limits of phase window (phase is 0->1)

        if autoroll, will roll the profile to center the pulse (makes some edge effects easier)

        """
        self.pfd = pfd(pfdfile)
        self.pfdfile = pfdfile
        self.pointing = Table.read(pointingfile)
        self.pointingfile = pointingfile
        self.RA = Longitude(
            np.interp(
                self.pfd.mid_secs, self.pointing["OFFS_SUB"], self.pointing["RA_SUB"]
            )
            * u.deg
        )
        self.Dec = Latitude(
            np.interp(
                self.pfd.mid_secs, self.pointing["OFFS_SUB"], self.pointing["DEC_SUB"]
            )
            * u.deg
        )
        # is it a RA scan or Dec scan
        if (self.RA.max() - self.RA.min()) * np.cos(self.Dec.mean()) > (
            self.Dec.max() - self.Dec.mean()
        ):
            self.scantype = "RA"
            self.offset = (self.RA - self.RA.mean()) * np.cos(self.Dec.mean())
        else:
            self.scantype = "Dec"
            self.offset = self.Dec - self.Dec.mean()

        self.freq = (self.pfd.hifreq + self.pfd.lofreq) / 2 * u.MHz
        if self.pfd.telescope == "GBT":
            self.diameter = 100 * u.m
        else:
            self.diameter = None

        if profilefile is not None:
            self.profile = PulseProfile(profilefile)
        else:
            self.profile = None

        # Axis 0 = subints, axis 1 = channels, axis 2 = phase bins
        self.pfd.dedisperse()
        self.subints = np.arange(self.pfd.profs.shape[0])
        self.channels = np.arange(self.pfd.profs.shape[1])
        self.phasebins = np.arange(self.pfd.profs.shape[2])
        self.phases = self.phasebins / (len(self.phasebins))
        self.rollvalue = 0

        # basic preparation
        # define off-pulse window
        self.set_offpulse(offpulse)
        # determine the bandpass
        self.get_bandpass()
        # remove the bandpass
        self.apply_bandpass()
        # remove the variable sub-integration power
        self.remove_subintpower()
        # roll the pulse to center it (if desired)
        if autoroll:
            self.autoroll()

        # fit a pulse profile if it is not already set
        if self.profile is None:
            self.profile = PulseProfile.fit_gaussian_profile(
                self.phases, np.nanmean(self.squeezed_profs, axis=0)
            )

    @property
    def telescope_fwhm(self):
        """
        OTF_Scan.telescope_fwhm

        property to return the telescope FWHM (based on observing frequency and diameter)
        """
        return (
            1.2 * (self.freq.to(u.m, equivalencies=u.spectral()) / self.diameter)
        ).to(u.deg, equivalencies=u.dimensionless_angles())

    @property
    def telescope_sigma(self):
        """
        OTF_Scan.telescope_sigma

        property to return the telescope sigma (based on observing frequency and diameter)
        """
        return self.telescope_fwhm / (2 * np.sqrt(2 * np.log(2)))

    def telescope_beam(self, x, offset, ampl):
        """
        beam = OTF_Scan.telescope_beam(x, offset, ampl)

        returns a beam profile evaluated at x with fixed FWHM and specified offset, amplitude
        """
        return gaussian(x, self.telescope_fwhm.to(u.deg).value, offset, ampl)

    def set_offpulse(self, offpulse):
        """
        OTF_Scan.set_offpulse(offpulse)
        offpulse will be pairs [min,max] where min and max are limits of phase window (phase is 0->1)
        turns it into an array of indices
        """
        if offpulse is not None and len(offpulse) > 0:
            good = np.zeros(len(self.phases), dtype=bool)
            for minphase, maxphase in offpulse:
                good = good | (
                    (self.phasebins >= (minphase * len(self.phases)))
                    & (self.phasebins <= (maxphase * len(self.phases)))
                )
            self.offpulse = np.where(good)[0]
        else:
            self.offpulse = self.phasebins

    def get_bandpass(self):
        """
        OTF_Scan.get_bandpass()
        determine the bandpass by taking the median of the off-pulse data averaging across time, phase
        """
        self.bandpass = np.nanmedian(
            self.pfd.profs[:, :, self.offpulse], axis=(0, 2), keepdims=True
        )

    def apply_bandpass(self):
        """
        OTF_Scan.apply_bandpass()
        subtract away the bandpass and take the mean across the frequency axis (squeezed)
        """

        self.squeezed_profs = np.nanmean(
            self.pfd.profs - self.bandpass, axis=1
        ).squeeze()

    def remove_subintpower(self):
        """
        OTF_Scan.remove_subintpower()
        calculate the power for each subintegration by taking the off-pulse median across phase
        and subtract it away
        """
        self.subint_power = np.nanmedian(
            self.squeezed_profs[:, self.offpulse], axis=(1), keepdims=True
        )
        self.squeezed_profs -= self.subint_power

    def roll(self):
        """
        OTF_Scan.roll()

        roll the phases (including offpulse and pulse profile object) to center the pulse
        uses the value of OTF_Scan.rollvalue
        """
        self.pfd.profs = np.roll(self.pfd.profs, self.rollvalue, axis=-1)
        self.squeezed_profs = np.roll(self.squeezed_profs, self.rollvalue, axis=-1)
        self.offpulse += self.rollvalue
        self.offpulse[self.offpulse < 0] += len(self.phasebins)
        if self.profile is not None:
            self.profile.pha += self.rollvalue / (len(self.phasebins))

    def autoroll(self):
        """
        OTF_Scan.autoroll()

        determine the roll to center the pulse based on the max of the profile (averaged across time)
        and apply it
        """
        profile = np.nanmedian(self.squeezed_profs, axis=0)
        self.rollvalue = int(
            len(self.phases) / 2 - self.phasebins[profile == profile.max()].mean()
        )
        self.roll()

    def iterate_flattening(self, update_profile=True, threshold=0.01):
        """
        OTF_Scan.iterate_flattening(update_profile=True, threshold=0.01)

        redo the various flattening steps (across bandpass, subintegrations)
        after recalculating off-pulse based on when the analytic profile is < threshold * max(analytic profile)

        if update_profile will recalculate the profile fit
        """
        self.offpulse = np.where(
            self.profile(self.phases) < threshold * self.profile(self.phases).max()
        )[0]
        self.get_bandpass()
        self.apply_bandpass()
        self.remove_subintpower()

        if update_profile:
            self.profile = PulseProfile.fit_gaussian_profile(
                self.phases, np.nanmean(self.squeezed_profs, axis=0)
            )

    def rolling_SNR(self, l=3):
        """
        SNR = OTF_Scan.rolling_SNR(l=3)

        computing a moving average SNR with a window that goes from -l to l subintegrations
        returns the SNR
        """

        SNR = np.zeros(len(self.subints))
        # SNR_err = np.zeros(len(self.subints))
        template = self.profile(self.phases, normalize=True)
        kernel = np.ones((2 * l + 1, 1))
        I = scipy.ndimage.convolve(self.squeezed_profs, kernel, mode="nearest")
        SNR = np.nansum(I * template, axis=1) / np.nansum(template ** 2)
        return SNR

    def fit_offset(self, l=3):
        """
        offset, ampltude = OTF_Scan.fit_offset(l=3)

        compute the rolling SNR and fit the telescope beam to it (with offset, amplitude free)
        return the offset, amplitude
        """
        SNR = self.rolling_SNR(l=l)
        good = ~np.isnan(SNR)  # & (SNR_err > 0)
        x = self.offset.deg[good]
        y = SNR[good]
        # dy = SNR_err[good]
        popt, pcov = curve_fit(
            self.telescope_beam,
            x,
            y,
            # sigma=dy,
            p0=[x[y == y.max()].mean(), y.max()],
            # absolute_sigma=True,
        )
        return popt[0] * u.deg, popt[1]  # , np.sqrt(pcov[0]) * u.deg, np.sqrt(pcov[1])


class OTF_ScanSet:
    def __init__(
        self,
        pfds,
        pointings,
        profilefiles=None,
        offpulses=None,
        autoroll=True,
        iterate=True,
    ):
        """
        OTF_ScanSet(pfds,
        pointings,
        profilefiles=None,
        offpulses=None,
        autoroll=True,
        iterate=True)

        load in a RA scan and Dec scan and do basic flattening

        pfds = [RA_pfd, Dec_pfd]
        pointings = [RA_pointing, Dec_pointing]
        profilefiles = [RA_profile, Dec_profile]
        offpulses = [RA_offpulse, Dec_offpulse]
        """

        if profilefiles is None or len(profilefiles)==0:
            profilefiles = [None, None]
        if offpulses is None or len(offpulses)==0:
            offpulses = [None, None]
        self.scan_RA = OTF_Scan(
            pfds[0],
            pointings[0],
            profilefile=profilefiles[0],
            offpulse=offpulses[0],
            autoroll=autoroll,
        )
        self.scan_Dec = OTF_Scan(
            pfds[1],
            pointings[1],
            profilefile=profilefiles[1],
            offpulse=offpulses[1],
            autoroll=autoroll,
        )
        assert self.scan_RA.scantype == "RA", "Scan '%s' has scantype='%s', not RA" % (
            pointings[0],
            self.scan_RA.scantype,
        )
        assert (
            self.scan_Dec.scantype == "Dec"
        ), "Scan '%s' has scantype='%s', not Dec" % (
            pointings[1],
            self.scan_Dec.scantype,
        )
        if iterate:
            if profilefiles[0] is None:
                self.scan_RA.iterate_flattening(update_profile=True)
            else:
                self.scan_RA.iterate_flattening(update_profile=False)
            if profilefiles[1] is None:
                self.scan_Dec.iterate_flattening(update_profile=True)
            else:
                self.scan_Dec.iterate_flattening(update_profile=False)

    def fit_offsets(self, l=5, plot=False):
        """
        position = OTF_ScanSet(l=5, plot=False)

        compute the best-fit position from the OTF Scans through RA and Dec
        if desired also make a plot
        """
        (
            self.offset_RA,
            self.ampl_RA,
            # self.offset_RA_err,
            # self.ampl_RA_err,
        ) = self.scan_RA.fit_offset(l=l)
        (
            self.offset_Dec,
            self.ampl_Dec,
            # self.offset_Dec_err,
            # self.ampl_Dec_err,
        ) = self.scan_Dec.fit_offset(l=l)
        position_RA = self.scan_RA.RA.mean() + self.offset_RA / np.cos(
            self.scan_RA.Dec.mean()
        )
        position_Dec = self.scan_Dec.Dec.mean() + self.offset_Dec
        if plot is not None and plot:
            SNR_RA = self.scan_RA.rolling_SNR(l=l)
            SNR_Dec = self.scan_Dec.rolling_SNR(l=l)
            plt.clf()
            plt.errorbar(self.scan_RA.offset.value, SNR_RA, fmt="o-", label="RA")
            plt.plot(
                self.scan_RA.offset,
                self.scan_RA.telescope_beam(
                    self.scan_RA.offset.value, self.offset_RA.value, self.ampl_RA
                ),
                "--",
                label="RA Beam",
            )
            plt.errorbar(
                self.scan_Dec.offset.value,
                SNR_Dec,
                # yerr=SNR_Dec_err,
                fmt="s-",
                label="Dec",
            )
            plt.plot(
                self.scan_Dec.offset,
                self.scan_Dec.telescope_beam(
                    self.scan_Dec.offset.value, self.offset_Dec.value, self.ampl_Dec
                ),
                "--",
                label="Dec Beam",
            )
            plt.xlabel("Position Offset (deg)")
            plt.ylabel("SNR")
            plt.legend()
            if isinstance(plot, str):
                plt.savefig(plot)
                print(f"Output saved to {plot}")

        return SkyCoord(position_RA, position_Dec)
