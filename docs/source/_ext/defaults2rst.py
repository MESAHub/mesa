import re

# regex to strip leading comment
STRIP = re.compile(r'^\s*\!\s?')

# files that should be stripped
DEFAULTS_FILES = (
    'reference/astero_search',
    'reference/astero_pgstar',
    'eos/defaults',
    'kap/defaults',
    'reference/binary_job',
    'reference/binary_controls',
    'reference/star_job',
    'reference/controls',
    'reference/pgstar',
    'reference/format',
)

def defaults2rst(app, docname, source):
    """
    Transform defaults files to ReST by stripping leading comment characters
    """
    if docname in DEFAULTS_FILES:
        src = source[0]
        lines = src.split('\n')
        rst_lines = [STRIP.sub('', line) for line in lines]
        source[0] = '\n'.join(rst_lines)
    return

def setup(app):
    app.connect("source-read", defaults2rst)
