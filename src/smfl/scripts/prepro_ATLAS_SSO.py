import sys
from pathlib import Path

import pandas as pd
from argparse import ArgumentParser


def split_by_night(df, jd_col="JD_lc", gap=0.2):
    """
    Split dataframe into nightly groups.
    """
    df = df.sort_values(jd_col).reset_index(drop=True)
    df["jd_gap"] = df[jd_col].diff()
    df["n_obs"] = (
        df["jd_gap"] > gap
    ).cumsum() + 1
    return df


def save_separated(
    df,
    stem,
    suffix,
    output_dir,
    ndatamin,
):

    if len(df) == 0:

        return {
            "initial": 0,
            "used": 0,
            "removed": 0,
        }

    initial_n = len(df)

    df = split_by_night(df)

    unique_nights = sorted(df["n_obs"].unique())

    counter = 1

    used_total = 0

    removed_total = 0

    for nobs in unique_nights:

        df_night = df[df["n_obs"] == nobs]

        ndata = len(df_night)

        # ======================================================
        # skip small datasets
        # ======================================================

        if ndata < ndatamin:

            removed_total += ndata

            print(
                f"skip {suffix} night {nobs}: "
                f"{ndata} points "
                f"(< {ndatamin})"
            )

            continue

        # ======================================================
        # save
        # ======================================================

        outfile = output_dir / (
            f"{stem}_{suffix}_{counter:03d}.csv"
        )

        df_night.to_csv(outfile, index=False)

        print(f"saved: {outfile}")

        used_total += ndata

        counter += 1

    return {
        "initial": initial_n,
        "used": used_total,
        "removed": removed_total,
    }


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "input",
        help="Input CSV file"
    )
    parser.add_argument(
        "--separate",
        action="store_true",
        help="Separate nightly data using JD_lc"
    )
    parser.add_argument(
        "--Ndatamin",
        type=int,
        default=3,
        help="Minimum number of data points per night"
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="Output directory"
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    # ==========================================================
    # output directory
    # ==========================================================
    
    if args.outdir is None:
    
        outdir = input_file.parent
    
    else:
        outdir = Path(args.outdir)
        outdir.mkdir(
            parents=True,
            exist_ok=True)

    # ==========================================================
    # read csv
    # ==========================================================

    df = pd.read_csv(input_file)

    # split by filter
    df_o = df[df["obs"].str.endswith("o")]

    df_c = df[df["obs"].str.endswith("c")]

    stem = input_file.stem

    # ==========================================================
    # normal save
    # ==========================================================

    if not args.separate:

        output_o = outdir / (
            f"{stem}_o.csv"
        )
        
        output_c = outdir / (
            f"{stem}_c.csv"
        )

        df_o.to_csv(output_o, index=False)

        df_c.to_csv(output_c, index=False)

        print(f"saved: {output_o}")

        print(f"saved: {output_c}")

        # summary
        print("\n================ SUMMARY ================")

        print(
            f"o band : initial = {len(df_o)}, "
            f"removed = 0"
        )

        print(
            f"c band : initial = {len(df_c)}, "
            f"removed = 0"
        )

    # ==========================================================
    # nightly separated save
    # ==========================================================

    else:

        stats_o = save_separated(
            df=df_o,
            stem=stem,
            suffix="o",
            output_dir=outdir,
            ndatamin=args.Ndatamin,
        )

        stats_c = save_separated(
            df=df_c,
            stem=stem,
            suffix="c",
            output_dir=outdir,
            ndatamin=args.Ndatamin,
        )

        # ======================================================
        # final summary
        # ======================================================

        print("\n================ SUMMARY ================")

        print(
            f"o band : "
            f"initial = {stats_o['initial']}, "
            f"used = {stats_o['used']}, "
            f"removed = {stats_o['removed']}"
        )

        print(
            f"c band : "
            f"initial = {stats_c['initial']}, "
            f"used = {stats_c['used']}, "
            f"removed = {stats_c['removed']}"
        )


if __name__ == "__main__":
    main()

