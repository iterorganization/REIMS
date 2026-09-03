from jinja2 import Environment
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('in_file', type=Path, help='input jinja2 file')
arg = parser.parse_args()
j_env = Environment(trim_blocks=True,lstrip_blocks=True)

with open(arg.in_file) as f:
    j_out = j_env.from_string(f.read()).render(place_holder='for future')
    with open(arg.in_file.with_suffix(''),'w') as out_f:
       out_f.write(j_out)
