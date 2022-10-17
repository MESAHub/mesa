from __future__ import print_function

import re
import sys


def search(filename, checks, summary=False):
    count = 0
    with open(str(filename), "r") as f:
        for ldx, line in enumerate(f):
            if any(line.startswith(i) for i in ["! ", "c ", "C ", "* "]):
                # Skip comment lines
                continue
            for c in checks:
                if "!" in line and not "!$" in line:
                    l = line[: line.index("!")]
                else:
                    l = line
                lcheck = l.lower()
                x = c(lcheck)
                if x is not None:
                    count += 1
                    if not summary:
                        print(line.strip())
                        print(str(filename), ":", str(ldx + 1), ":", x)
                        print()
    if summary and count:
        print(str(filename), "Count:", count)


def check_float(line):
    if "float(" in line:
        return "Don't use float() use dble()"
    return None


def check_pow(line):
    # Look for 3**5
    if re.search("[a-zA-Z0-9\)]\*\*[a-zA-Z0-9]", line):
        return "Found ** use, use powX() instead"
    return None


def check_real_op(line):
    # Lots of code has 1.+2.
    checks = ["[0-9]\.\*", "[0-9]\.\+", "[0-9]\.\-", "[0-9]\.\/"]
    found = []
    for c in checks:
        if re.search(c, line):
            found.append(c)
    if len(found):
        return "Single precision number found "
    return None


def check_real_exp(line):
    # Look for 1e+1, 1e-1, 1e1
    if (
            re.search("[0-9][eE][+]?[-]?[0-9]", line)
            and "write" not in line
            and "format(" not in line
    ):
        return "Found use of exponent E, use D instead"
    return None


def check_real_d(line):
    # Look for 1.5 but not 1.5d0
    if "write" not in line and "format(" not in line:
        for i in re.split(
                " |\+|\-|\=|\*|\/", line
        ):  # Split up string into things approximately like a number
            if len(i) and re.search(
                    "\d+[.](?!.*[_Dd]).*", i
            ):  # test if missing double precision qualifier
                return "Missing D on float"
    return None


def check_real(line):
    # Look for declaring things real and not real(dp)
    if "real " in line or "real," in line:
        return "Declared real use real(dp) instead"
    return None


def check_dp(line):
    if "double precision" in line:
        return "Found double precision use real(dp) instead"
    return None


allchecks = [
    check_float,
    check_pow,
    check_real_op,
    check_real_exp,
    check_real_d,
    check_dp,
]

if __name__ == "__main__":
    files = sys.argv[1:]
    s = False
    if "-s" in files:
        s = True
        files = [i for i in files if not ("-s" in i)]
    for f in files:
        search(f, allchecks, s)
