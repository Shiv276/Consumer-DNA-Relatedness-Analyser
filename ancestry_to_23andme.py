import re
from pathlib import Path


def process(filename, output_prefix):
    output_filename = Path(f"{output_prefix}.txt")

    ancestry_pattern = r"^(rs\d+)\s+(\d+|X|Y|MT)\s+(\d+)\s+([ATGC])\s+([ATGC])"

    already_pattern = r"^(rs\d+)\s+(\d+|X|Y|MT)\s+(\d+)\s+([ATGC]{2})"

    with open(filename, "r") as f1:
        with open(output_filename, "w") as w2:

            for line in f1:

                if line.startswith("#"):
                    continue

                elif line.startswith("rsid\tchromosome\tposition\tallele1\tallele2"):
                    continue

                elif re.search(ancestry_pattern, line):
                    line = re.sub(ancestry_pattern, r"\1\t\2\t\3\t\4\5", line)
                    w2.write(line)

                elif re.search(already_pattern, line):
                    w2.write(line)

    return output_filename