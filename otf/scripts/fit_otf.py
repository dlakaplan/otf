#!/usr/bin/env python
import sys
import os
import argparse
import otf


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--pfd", type=str, nargs=2, required=True, help="Input PFD files"
    )
    parser.add_argument(
        "--pointing", type=str, nargs=2, required=True, help="Input pointing tables",
    )
    parser.add_argument(
        "--profile",
        type=str,
        nargs=2,
        help="Input pulse profiles (output of pygaussfit)",
    )
    parser.add_argument(
        "-l",
        "--length",
        default=5,
        type=int,
        help="Length of sliding window (subintegrations)",
    )
    parser.add_argument(
        "--noroll", action="store_false", help="Do not auto-roll the pulses"
    )

    parser.add_argument(
        "--offpulse",
        default=None,
        nargs=2,
        help="Off-pulse windows (e.g., '[[0,0.2],[0.8,1.0]]' for each scan)",
    )
    parser.add_argument("--plot", default=None, help="Name of output plot")

    args = parser.parse_args()
    offpulses = list(map(eval, args.offpulse))

    o = otf.OTF_ScanSet(
        args.pfd,
        args.pointing,
        profilefiles=args.profile,
        offpulses=list(map(eval, args.offpulse)),
        autoroll=~args.noroll,
    )

    p = o.fit_offsets(l=args.length, plot=args.plot)
    print(
        "From {} and {} determined best-fit position of {}".format(
            args.pfd[0], args.pfd[1], p.to_string("hmsdms")
        )
    )


if __name__ == "__main__":
    main()
