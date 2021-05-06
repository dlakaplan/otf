#!/usr/bin/env python
import sys
import os
import argparse
import otf


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("files", type=str, nargs="+", help="Input PSRFITS file(s)")
    parser.add_arcument(
        "-s",
        "--suffix",
        type=str,
        default="pointing",
        help="Suffix to append for output files",
    )
    parser.add_argument(
        "--format", choices=["ecsv", "fits"], help="Format of output table"
    )

    args = parser.parse_args()

    for filename in args.files:
        outputfile = otf.extract_pointing(
            filename, suffix=args.suffix, format=args.format
        )
        print(f"Extracted pointing information from {filename} to {outputfile}")


if __name__ == "__main__":
    main()
