import glob
import os

# Checks whether a test case has a owner and warns for test cases with only
# one owner

MESA_DIR = os.environ["MESA_DIR"]

CODEOWNERS = os.path.join(MESA_DIR, "CODEOWNERS")

STAR = os.path.join(MESA_DIR, "star", "test_suite")
BINARY = os.path.join(MESA_DIR, "binary", "test_suite")
ASTERO = os.path.join(MESA_DIR, "astero", "test_suite")


def get_test_cases(folder):
    # Get test cases
    x = list(glob.glob(os.path.join(folder, "*/")))
    return set(
        [i.removeprefix(MESA_DIR).removesuffix("/").removeprefix("/") for i in
         x]
    )


def get_do1(folder):
    test_cases = []
    with open(os.path.join(folder, "do1_test_source"), "r") as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith("#") or not len(line) or line == "return":
                continue
            test_cases.append(line.split()[1])
    
    return set(test_cases)


def parse_codeowners(filename):
    with open(filename) as f:
        lines = f.readlines()
    
    result = {
        "star": {},
        "binary": {},
        "astero": {},
        "module": {},
    }
    for line in lines:
        if line.startswith("#") or len(line.strip()) == 0:
            continue
        
        case = line.split()
        
        if case[0].startswith("star/test_suite"):
            result["star"][case[0]] = case[1:]
        elif case[0].startswith("binary/test_suite"):
            result["binary"][case[0]] = case[1:]
        elif case[0].startswith("astero/test_suite"):
            result["astero"][case[0]] = case[1:]
        else:
            result["module"][case[0]] = case[1:]
    
    return result


def print_section(header):
    """Display output section header"""
    print(f"\n\n*** {header} ***\n")


def print_options(options):
    """Print a set of options"""
    for o in sorted(options):
        print(f"   {o}")


# Things with not in CODEOWNERS file
def not_listed(cases, code):
    listed = set(code.keys())
    result = cases - listed
    if len(result):
        print_section("Test cases not in CODEOWNERS")
        print_options(result)


def no_longer_exists(cases, code):
    listed = set(code.keys())
    result = listed - cases
    if len(result):
        print_section("Test cases in CODEOWNERS but not do1_test_source")
        print_options(result)


# Check number of owners
def check_owners(code, num=0):
    result = []
    for key, value in code.items():
        if len(value) == num:
            result.append(key)
    
    if len(result):
        print_section(f"Test cases with only {num} CODEOWNERS")
        print_options(result)


if __name__ == "__main__":
    star_cases = get_test_cases(STAR)
    binary_cases = get_test_cases(BINARY)
    astero_cases = get_test_cases(ASTERO)
    
    star_do1 = get_test_cases(STAR)
    binary_do1 = get_test_cases(BINARY)
    astero_do1 = get_test_cases(ASTERO)
    
    codeowner = parse_codeowners(CODEOWNERS)
    
    not_listed(star_cases, codeowner["star"])
    not_listed(binary_cases, codeowner["binary"])
    not_listed(astero_cases, codeowner["astero"])
    
    no_longer_exists(star_do1, codeowner["star"])
    no_longer_exists(binary_do1, codeowner["binary"])
    no_longer_exists(astero_do1, codeowner["astero"])
    
    check_owners(codeowner["star"], 0)
    check_owners(codeowner["binary"], 0)
    check_owners(codeowner["astero"], 0)
    
    check_owners(codeowner["star"], 1)
    check_owners(codeowner["binary"], 1)
    check_owners(codeowner["astero"], 1)
