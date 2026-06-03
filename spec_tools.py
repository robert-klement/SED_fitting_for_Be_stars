import datetime as _dt

import numpy as np
import astropy.io.fits as fits
import astropy.units as u
from astropy.time import Time
from astropy.modeling import fitting as ap_fitting


MJD_JD2000 = 51544.5


def _to_quantity(values, default_unit):
    if isinstance(values, u.Quantity):
        return values
    return np.asarray(values) * default_unit


class Spectrum1D:
    """Minimal 1D spectrum container used by scripts in this repository."""

    def __init__(self, spectral_axis, flux):
        self.spectral_axis = _to_quantity(spectral_axis, u.AA)
        self.flux = _to_quantity(flux, u.one)
        if self.spectral_axis.shape != self.flux.shape:
            raise ValueError(
                "spectral_axis and flux must have the same shape, "
                f"got {self.spectral_axis.shape} and {self.flux.shape}"
            )

    @property
    def wavelength(self):
        return self.spectral_axis

    def __len__(self):
        return len(self.flux)

    def __truediv__(self, other):
        if isinstance(other, Spectrum1D):
            denom = other.flux
        else:
            denom = other
        return Spectrum1D(spectral_axis=self.spectral_axis, flux=self.flux / denom)


class SpectralRegion:
    """Simple spectral region with lower/upper bounds."""

    def __init__(self, lower, upper):
        self.lower = _to_quantity(lower, u.AA)
        self.upper = _to_quantity(upper, u.AA)
        if self.lower > self.upper:
            self.lower, self.upper = self.upper, self.lower


def _normalize_regions(regions):
    if regions is None:
        return []
    if isinstance(regions, SpectralRegion):
        return [regions]
    if (
        isinstance(regions, tuple)
        and len(regions) == 2
        and not isinstance(regions[0], SpectralRegion)
    ):
        return [SpectralRegion(regions[0], regions[1])]

    out = []
    for item in regions:
        if isinstance(item, SpectralRegion):
            out.append(item)
        elif isinstance(item, tuple) and len(item) == 2:
            out.append(SpectralRegion(item[0], item[1]))
    return out


def _region_mask(axis, regions):
    regs = _normalize_regions(regions)
    if not regs:
        return np.ones(axis.shape, dtype=bool)

    mask = np.zeros(axis.shape, dtype=bool)
    for reg in regs:
        mask |= (axis >= reg.lower) & (axis <= reg.upper)
    return mask


def _flux_values(spec_or_flux):
    flux = spec_or_flux.flux if hasattr(spec_or_flux, "flux") else spec_or_flux
    if hasattr(flux, "value"):
        return np.asarray(flux.value, dtype=float)
    return np.asarray(flux, dtype=float)


def snr_derived(spec_or_flux):
    """SNR estimate from median signal and point-to-point noise."""
    flux = _flux_values(spec_or_flux)
    flux = flux[np.isfinite(flux)]
    if flux.size < 3:
        return np.nan

    noise = np.nanstd(np.diff(flux)) / np.sqrt(2.0)
    if not np.isfinite(noise) or noise <= 0:
        return np.inf

    signal = np.nanmedian(flux)
    return float(signal / noise)


class _LinearContinuum:
    def __init__(self, slope, intercept, unit):
        self.slope = float(slope)
        self.intercept = float(intercept)
        self.unit = unit

    def __call__(self, x):
        if isinstance(x, u.Quantity):
            xv = x.to_value(u.AA)
        else:
            xv = np.asarray(x, dtype=float)
        return (self.slope * xv + self.intercept) * self.unit


def fit_continuum(spectrum, window=None, exclude_regions=None, fitter=None):
    """Fit a simple linear continuum, compatible with old specutils calls."""
    del fitter  # kept for call-site compatibility
    axis = spectrum.spectral_axis
    x = np.asarray(axis.to_value(u.AA), dtype=float)
    y = _flux_values(spectrum)

    good = np.isfinite(x) & np.isfinite(y)
    if window is not None:
        good &= _region_mask(axis, window)
    if exclude_regions is not None:
        good &= ~_region_mask(axis, exclude_regions)

    if np.count_nonzero(good) < 2:
        slope = 0.0
        intercept = float(np.nanmedian(y[np.isfinite(y)])) if np.any(np.isfinite(y)) else 1.0
    else:
        slope, intercept = np.polyfit(x[good], y[good], 1)

    return _LinearContinuum(slope, intercept, spectrum.flux.unit)


def fit_lines(spectrum, model, fitter=None, window=None):
    """Fit a model/compound model to a 1D spectrum."""
    if fitter is None:
        fitter = ap_fitting.LevMarLSQFitter()

    axis = spectrum.spectral_axis
    flux = spectrum.flux
    good = np.isfinite(axis.value) & np.isfinite(flux.value)
    if window is not None:
        good &= _region_mask(axis, window)

    if np.count_nonzero(good) < 3:
        raise ValueError("Not enough valid points for line fit.")

    return fitter(model, axis[good], flux[good])


def estimate_line_parameters(spectrum, model):
    """Compatibility no-op for legacy scripts."""
    del spectrum
    return model


def equivalent_width(spectrum, regions=None):
    """Compute equivalent width in Angstrom."""
    axis = spectrum.spectral_axis
    wav = np.asarray(axis.to_value(u.AA), dtype=float)
    flux = _flux_values(spectrum)

    regs = _normalize_regions(regions)
    if not regs:
        regs = [None]

    ew_total = 0.0
    for reg in regs:
        if reg is None:
            mask = np.ones(wav.shape, dtype=bool)
        else:
            mask = (axis >= reg.lower) & (axis <= reg.upper)
        good = mask & np.isfinite(wav) & np.isfinite(flux)
        if np.count_nonzero(good) < 2:
            continue
        ew_total += np.trapezoid(1.0 - flux[good], wav[good])

    return ew_total * u.AA


def _build_wavelength_axis(header, npts):
    crval1 = header.get("CRVAL1")
    cdelt1 = header.get("CDELT1", header.get("CD1_1"))
    crpix1 = float(header.get("CRPIX1", 1.0))

    if crval1 is None or cdelt1 is None:
        raise KeyError("Missing CRVAL1/CDELT1 (or CD1_1) in FITS header.")

    pix = np.arange(npts, dtype=float) + 1.0
    return crval1 + (pix - crpix1) * cdelt1


def _parse_date_obs_to_mjd(date_obs):
    if not date_obs:
        return None

    text = str(date_obs).strip()
    if not text:
        return None

    for fmt in ("isot", "iso"):
        try:
            return float(Time(text, format=fmt, scale="utc").mjd)
        except Exception:
            pass

    if "T" in text:
        dpart, tpart = text.split("T", 1)
    else:
        dpart, tpart = text, "00:00:00"

    try:
        if "-" in dpart:
            y, m, d = dpart.split("-")
            dt_obj = _dt.datetime(
                int(y), int(m), int(d),
                int(tpart.split(":")[0]),
                int(tpart.split(":")[1]),
                int(float(tpart.split(":")[2])),
            )
            return float(Time(dt_obj, scale="utc").mjd)
    except Exception:
        pass

    try:
        if "/" in dpart:
            d, m, y = dpart.split("/")
            dt_obj = _dt.datetime(
                int(y), int(m), int(d),
                int(tpart.split(":")[0]),
                int(tpart.split(":")[1]),
                int(float(tpart.split(":")[2])),
            )
            return float(Time(dt_obj, scale="utc").mjd)
    except Exception:
        pass

    return None


def _extract_mjd_from_header(header):
    if "MJD-OBS" in header:
        return float(header["MJD-OBS"])
    if "MJD" in header:
        return float(header["MJD"])
    if "JD" in header:
        return float(header["JD"]) - 2400000.5
    if "HJD" in header:
        # Some files store HJD directly; retain historical behavior.
        return float(header["HJD"])

    for key in ("DATE-OBS", "FRAME"):
        if key in header:
            mjd = _parse_date_obs_to_mjd(header[key])
            if mjd is not None:
                return mjd

    return MJD_JD2000


def loadfits(fitsfile, data_header=0, ARCES_reduced_by_Peter=False):
    del ARCES_reduced_by_Peter  # compatibility argument retained

    with fits.open(fitsfile) as hdul:
        data = np.asarray(hdul[data_header].data).squeeze()
        if data.ndim != 1:
            raise ValueError(f"Expected 1D spectrum in {fitsfile}, got shape {data.shape}.")

        header = hdul[0].header
        flux = data
        wl = _build_wavelength_axis(header, flux.size)

        if "WLSHIFT" in header:
            wl = wl + float(header["WLSHIFT"])

        MJD = _extract_mjd_from_header(header)
        dateobs = header.get("DATE-OBS", header.get("FRAME", ""))
        datereduc = header.get("IRAF-TLM", header.get("DATE", ""))

    return wl, flux, MJD, dateobs, datereduc, fitsfile


def loadfits_rivi(fitsfile):
    with fits.open(fitsfile) as hdul:
        data = np.asarray(hdul[0].data).squeeze()
        if data.ndim != 1:
            raise ValueError(f"Expected 1D spectrum in {fitsfile}, got shape {data.shape}.")

        header = hdul[0].header
        flux = data
        wl = _build_wavelength_axis(header, flux.size)

        if "WLSHIFT" in header:
            wl = wl + float(header["WLSHIFT"])

        MJD = None
        history = header.get("HISTORY", [])
        if isinstance(history, str):
            history = [history]

        for i, line in enumerate(history):
            if "DTIME" in str(line) and i + 2 < len(history):
                try:
                    MJD = float(str(history[i + 2])[2:26])
                except Exception:
                    MJD = None
                break

        if MJD is None:
            MJD = _extract_mjd_from_header(header)

    return wl, flux, MJD


def get_spec_snr(wav, flux, region_limits):
    mask = (wav >= region_limits[0]) & (wav <= region_limits[1])
    if np.count_nonzero(mask) < 3:
        return np.nan

    spec = type(
        "_SimpleSpec",
        (),
        {
            "spectral_axis": wav[mask] * u.AA,
            "flux": flux[mask] * u.one,
        },
    )()
    return snr_derived(spec)
