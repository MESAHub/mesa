import os
import glob
from collections import defaultdict

import check_test_suite_owners as cc

# Lists test cases for each name


codeowner = cc.parse_codeowners(cc.CODEOWNERS)

def list_by_author(codeowner):

    result = defaultdict(list)

    for module in ['star','binary','astero']:
        code = codeowner[module]
        for key, value in  code.items():
            if len(value) == 0:
                continue

            for name in value:
                result[name].append(key)

        
    for key,value in result.items():
        cc.print_section(f'Number of test cases owned by {key} = {len(value)}')
        cc.print_options(value)


if __name__ == "__main__":
    list_by_author(codeowner)
