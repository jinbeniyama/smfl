import numpy as np
import pandas as pd
from argparse import ArgumentParser as ap
from astroquery.jplhorizons import Horizons
from astropy import units as u
import warnings

warnings.filterwarnings("ignore")

C_SPEED = 3e5  # km/s
AU2KM = (1 * u.au).to("km").value
SEC2DAY = 24.0 * 3600.


def calc_ltday_vectorized(obj, code, t_jds):
    """Calculate light travel time [day] for an array of JDs."""
    try:
        # Querying JPL Horizons with a list/array of epochs is much faster
        ast = Horizons(id=obj, location=code, epochs=list(t_jds))
        eph = ast.ephemerides()

        # Extract delta values (ensuring alignment with the epochs)
        delta = eph["delta"].values  # AU
        lt_day = (delta * AU2KM) / C_SPEED / SEC2DAY

        return lt_day

    except Exception as e:
        print(f"Error fetching Horizons data for the requested epochs: {e}")
        # Fallback to single-item processing if vectorization fails due to formatting
        print("Attempting row-by-row fallback calculation...")
        return None


def calc_ltday_single(obj, code, t_jd):
    """Fallback function for a single JD calculation."""
    try:
        ast = Horizons(id=obj, location=code, epochs=t_jd)
        eph = ast.ephemerides()
        delta = eph["delta"][0]
        return (delta * AU2KM) / C_SPEED / SEC2DAY
    except Exception as e:
        print(f"Error fetching Horizons data for JD {t_jd}: {e}")
        return np.nan


def get_args():
    parser = ap(description="Light-time correction for parsed Lightcurve CSV")
    parser.add_argument(
        "obj",
        type=str,
        help="Object name/number for JPL Horizons"
    )
    parser.add_argument(
        "lc",
        type=str,
        help="Input CSV file"
    )
    parser.add_argument(
        "out",
        type=str,
        help="Output filename"
    )
    parser.add_argument(
        "--code",
        type=str,
        default="500",
        help="MPC code (default: 500)"
    )
    parser.add_argument(
        "--f_max",
        type=float,
        default=None,
        help="Maximum flux to save"
    )
    parser.add_argument(
        "--snr_min",
        type=float,
        default=None,
        help="Minimum SNR to save"
    )
    
    # Optional column naming arguments for higher versatility
    parser.add_argument(
        "--col_jd",
        type=str,
        default="jd",
        help="Column name for Julian Date (default: jd)"
    )
    parser.add_argument(
        "--col_mag",
        type=str,
        default="mag",
        help="Column name for Magnitude (default: mag)"
    )
    parser.add_argument(
        "--col_magerr",
        type=str,
        default="magerr",
        help="Column name for Magnitude Error (default: magerr)"
    )
    parser.add_argument(
        "--col_band",
        type=str,
        default="band",
        help="Column name for Filter/Band (default: band)"
    )

    return parser.parse_args()


def main(args=None):

    if args is None:
        args = get_args()

    obj = args.obj
    code = args.code

    print(f"Loading input file: {args.lc}")
    df = pd.read_csv(args.lc)

    # Dynamically validate columns based on arguments
    required_cols = [args.col_jd, args.col_mag, args.col_magerr, args.col_band]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(
                f"Required column '{col}' missing from input CSV. "
                f"Available columns: {list(df.columns)}"
            )

    # Sort based on custom JD column
    df = df.sort_values(args.col_jd).reset_index(drop=True)

    # Session split
    df["session_gap"] = df[args.col_jd].diff()
    df["n_obs"] = (df["session_gap"] > 0.2).cumsum() + 1
    df["sdate"] = df.groupby("n_obs")[args.col_jd].transform(
        lambda x: str(int(np.min(x)))
    )

    # ==========================================================
    # Light-time correction
    # ==========================================================
    print("\nCalculating light-time correction...")
    
    # Try efficient vectorized batch request first
    lt_days = calc_ltday_vectorized(obj, code, df[args.col_jd].values)
    
    # If vectorized query fails, fall back to row-by-row iteration
    if lt_days is None:
        lt_days = []
        for i, row in df.iterrows():
            t_jd = row[args.col_jd]
            print(f"[{i+1}/{len(df)}] JD = {t_jd:.8f}")
            lt_day = calc_ltday_single(obj=obj, code=code, t_jd=t_jd)
            lt_days.append(lt_day)

    df["lt_day"] = lt_days
    df["t_jd_ltcor"] = df[args.col_jd] - df["lt_day"]

    if "session_gap" in df.columns:
        df = df.drop(columns=["session_gap"])

    # ==========================================================
    # Save combined file
    # ==========================================================
    out_all = args.out
    df.to_csv(out_all, index=False)
    print(f"\nSaved combined data to: {out_all}")

    # ==========================================================
    # Save session-separated files
    # ==========================================================
    unique_sessions = sorted(df["n_obs"].unique())

    if len(unique_sessions) > 1:
        for nobs in unique_sessions:
            df_temp = df[df["n_obs"] == nobs]
            fname = df_temp["sdate"].iloc[0]

            # Optional filtering checks
            if args.f_max and "flux_cor" in df_temp.columns:
                df_temp = df_temp[df_temp["flux_cor"] < args.f_max]

            if args.snr_min and "flux_cor" in df_temp.columns:
                df_temp = df_temp[
                    (df_temp["flux_cor"] / df_temp["fluxerr_cor"]) > args.snr_min
                ]

            out_sep = f"{args.out}_{fname}"
            df_temp.to_csv(out_sep, index=False)
            print(f"Saved session file: {out_sep}")


if __name__ == "__main__":
    main()
